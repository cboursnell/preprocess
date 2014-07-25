#!/usr/bin/ruby

## # # # # # # # # # # # # #
## mrna pipeline
##
## trimmomatic
## bayeshammer (optional)
## khmer
##
## created: 2014-05-27 Chris Boursnell (cmb211@cam.ac.uk)
##
## # # # # # # # # # # # # #

require 'rubygems'
#require 'preprocessor'
require 'set'
require 'cmd'
require 'which'
require 'bindeps'
include Which

class MalformedInputError < StandardError
end

class Preprocessor

  attr_accessor :input, :output, :phred
  attr_reader :data

  def initialize(output, verbose, threads=1, memory=4)
    @verbose = verbose
    @trim_jar = "bin/trimmomatic-0.32.jar"
    @khmer = which("normalize-by-median.py").first
    @hammer_path = "bin/SPAdes-3.1.0-Linux/bin/spades.py"
    @output_dir = output ? File.expand_path(output) : Dir.pwd
    @memory = memory
    @threads = threads
    @data = []
  end

  def load_input(input)
    if File.exist?(input)
      File.open("#{input}").each_line do |line|
        cols = line.chomp.split(",")
        if cols.size != 5
          raise MalformedInputError.new("Input file does not contain 5 columns")
        end
        raise RuntimeError.new("#{cols[0]} not found") if !File.exist?(cols[0])
        if cols[4].to_i != 1 and cols[4].to_i != 2
          raise RuntimeError.new("Pair should be 1 or 2")
        end
        @data << { :name => cols[0],
                   :file => File.expand_path(cols[1]),
                   :rep => cols[2].to_i,
                   :type => cols[3],
                   :pair => cols[4].to_i,
                   :current => File.expand_path(cols[1]) }
      end
      @paired = @data.reduce(0) {|max,v| max=[max,v[:pair]].max}
    else
      raise RuntimeError, "#{input} does not exist"
    end
  end

  def load_reads(left, right, name)
    @data = []
    rep = 1
    left.split(",").zip(right.split(",")).each do |a, b|
      # left
      @data << { :name => name,
                 :file => File.expand_path(a),
                 :rep => rep,
                 :type => name,
                 :pair => 1,
                 :current => File.expand_path(a) }
      # right
      @data << { :name => name,
                 :file => File.expand_path(b),
                 :rep => rep,
                 :type => name,
                 :pair => 2,
                 :current => File.expand_path(b) }
      rep += 1
    end
    @paired = 2
  end

  def gunzip
    @data.each do |info|
      if info[:file]=~/\.gz$/
        output_filename = File.basename(info[:file].split(".gz").first)
        output_filename = File.join(@output_dir, output_filename)
        File.open("#{output_filename}", "wb") do |out|
          Zlib::GzipReader.open(info[:file]) do |gz|
            out.write(gz.read)
          end
        end
        info[:current] = output_filename
      end
    end
  end

  def set_output(output_dir)
    @output_dir = File.expand_path(output_dir)
  end

  def detect_phred
    file = File.open(@data[0][:current])
    c = 0
    scores={}
    while c < 400
      line = file.readline
      if c % 4 == 3
        line.chars.each do |i|
          ascii_key = i.ord
          scores[ascii_key] ||= 0
          scores[ascii_key] += 1
        end
      end
      c += 1
    end
    max = scores.keys.max
    min = scores.keys.min
    @phred = -1
    if max <= 74
      @phred = 33
    elsif max <= 105
      @phred = 64
    else
      raise RuntimeError.new("Couldn't determine phred formatting")
    end
    file.close # is this correct?
    return @phred
  end

  def trimmomatic(minlen=40, windowsize=4, quality=15, trailing=15,
            leading=15, mismatches=2)
    if @paired==1
      @data.each_with_index do |a, i|
        puts "a = #{a}"
        puts "i = #{i}"
      end
    elsif @paired==2
      @data.each_with_index.each_slice(2) do |(a,i), (b,j)|
        outfile_left = "#{@output_dir}/#{a[:type]}_#{a[:rep]}-#{a[:pair]}.t.fq"
        outfile_right = "#{@output_dir}/#{b[:type]}_#{b[:rep]}-#{b[:pair]}.t.fq"
        outfileU_left = "#{@output_dir}/#{a[:type]}_#{a[:rep]}-#{a[:pair]}.tU.fq"
        outfileU_right = "#{@output_dir}/#{b[:type]}_#{b[:rep]}-#{b[:pair]}.tU.fq"
        trim_cmd = "java -jar #{@trim_jar} PE "
        trim_cmd << " -phred#{self.detect_phred} "
        trim_cmd << " -threads #{@threads} "
        trim_cmd << " #{a[:current]} #{b[:current]} "
        trim_cmd << " #{outfile_left} #{outfileU_left} "
        trim_cmd << " #{outfile_right} #{outfileU_right} "
        trim_cmd << " LEADING:#{leading} TRAILING:#{trailing} "
        trim_cmd << " SLIDINGWINDOW:#{windowsize}:#{quality} MINLEN:#{minlen}"
        @data[i][:current] = outfile_left
        @data[j][:current] = outfile_right
        @data[i][:unpaired] = outfileU_left
        @data[j][:unpaired] = outfileU_right
        if !File.exist?("#{outfile_left}")
          # puts trim_cmd if @verbose
          cmd = Cmd.new(trim_cmd)
          cmd.run
        else
          puts "trimmomatic already run on #{a[:file]}" if @verbose
        end
      end
    end
  end

  def construct_hammer_input
    if @paired==2
      left=[]
      right=[]
      single=[]
      @data.each_with_index.each_slice(2) do |(a,i), (b,j)|
        # build a yaml file
        # run spades.py on all files
        left << a[:current]
        right << b[:current]
        single << a[:unpaired] if a[:unpaired]
        single << b[:unpaired] if b[:unpaired]
      end
      yaml = "[\n  {\n    orientation: \"fr\",\n    type: \"paired-end\",\n"
      yaml << "    left reads: [\n"
      left.each_with_index do |left_read, j|
        yaml << "      \"#{left_read}\""
        yaml << "," if j < left.length-1
        yaml << "\n"
      end
      yaml << "    ],\n    right reads: [\n"
      right.each_with_index do |right_read, j|
        yaml << "      \"#{right_read}\""
        yaml << "," if j < right.length-1
        yaml << "\n"
      end
      yaml << "    ],\n  }"
      if single.size > 0
        yaml << ",\n  {\n"
        yaml << "    type: \"single\",\n    single reads: [\n"
        single.each_with_index do |single_read, j|
          yaml << "      \"#{single_read}\""
          yaml << "," if j < single.length-1
          yaml << "\n"
        end
        yaml << "    ]\n  }"
      end
      yaml << "\n]\n"
    end
    # puts yaml
    File.open("#{@output_dir}/dataset.yaml","w") {|io| io.write yaml}
  end

  def hammer
    # wget http://spades.bioinf.spbau.ru/release3.0.0/SPAdes-3.0.0-Linux.tar.gz
    # tar -xzf SPAdes-3.0.0-Linux.tar.gz
    # cd SPAdes-3.0.0-Linux/bin/
    # --only-error-correction runs only read error correction (no assembly)
    # --disable-gzip-output forces error correction not to compress the
    #                       corrected reads
    if !File.exist?(@hammer_path)
      wget_cmd = "wget http://spades.bioinf.spbau.ru/release3.1.0/"
      wget_cmd << "SPAdes-3.1.0-Linux.tar.gz -O bin/spades.tar.gz"
      tar_cmd = "tar xzf bin/spades.tar.gz -C bin"
      # wget
      download = Cmd.new(wget_cmd)
      download.run
      puts download.cmd if @verbose

      # tar
      unpack = Cmd.new(tar_cmd)
      unpack.run
      puts unpack.cmd if @verbose
    end
    construct_hammer_input

    cmd = "python #{@hammer_path} --dataset #{@output_dir}/dataset.yaml "
    cmd << " --only-error-correction "
    cmd << " --disable-gzip-output -m #{@memory} -t #{@threads}"
    cmd << " -o #{@output_dir}/hammer"

    if !Dir.exist?("#{@output_dir}/hammer/corrected")
      c = Cmd.new(cmd)
      c.run
      c.cmd if @verbose
    end
    section=nil
    corrected_left=[]
    corrected_right=[]
    corrected_single=[]
    File.open("#{@output_dir}/hammer/corrected/corrected.yaml").each_line do |line|
      if line =~ /left.reads/
        section = :left
      elsif line =~ /right.reads/
        section = :right
      elsif line =~ /single.reads/
        section = :single
      elsif line =~ /fastq/
        if section == :left
          filename = line.split(/\s+/).last
          corrected_left << filename
        elsif section == :right
          filename = line.split(/\s+/).last
          corrected_right << filename
        elsif section == :single
          filename = line.split(/\s+/).last
          corrected_single << filename
        end
      end
    end

    # match the output back with the original file
    @data.each_with_index.each_slice(2) do |(a,i), (b,j)|
      left = a[:current]
      right = b[:current]
      leftU = a[:unpaired] if a[:unpaired]
      rightU = b[:unpaired] if b[:unpaired]
      left = File.basename(left).split(".")[0..-2].join(".")
      right = File.basename(right).split(".")[0..-2].join(".")
      leftU = File.basename(leftU).split(".")[0..-2].join(".") if a[:unpaired]
      rightU = File.basename(rightU).split(".")[0..-2].join(".") if b[:unpaired]
      corrected_left.each do |corr|
        if corr =~ /#{left}/
          @data[i][:current] = corr
        end
      end
      corrected_right.each do |corr|
        if corr =~ /#{right}/
          @data[j][:current] = corr
        end
      end
      if corrected_single.size > 0
        corrected_single.each do |corr|
          if corr =~ /#{leftU}/
            @data[i][:unpaired] = corr
          end
          if corr =~ /#{rightU}/
            @data[j][:unpaired] = corr
          end
        end
      end
    end
    # loads output file location into @data[:current]
  end

  def khmer(kmer=23, cutoff=20, buckets=4)
    if @khmer.nil?
      raise RuntimeError.new("khmer normalize-by-median not installed")
    end
    x = (@memory/buckets*1e9).to_i
    # interleave the input files if paired
    pair = ""
    if @paired == 2
      @data.each_with_index.each_slice(2) do |(a,i), (b,j)|
        input_left = a[:current]
        input_right = b[:current]
        if input_left and input_right
          outfile = "#{@output_dir}/"
          outfile << "#{a[:type]}_#{a[:rep]}-#{a[:pair]}"
          outfile << ".in.fq"

          cmd = "paste #{input_left} #{input_right} | paste - - - - | "
          cmd << "awk -v FS=\"\\t\" -v OFS=\"\\n\" \'{print(\"@read\"NR\":1\","
          cmd << "$3,$5,$7,\"@read\"NR\":2\",$4,$6,$8)}\' > "
          cmd << " #{outfile}"
          if !File.exist?(outfile)
            c = Cmd.new(cmd)
            c.run
            puts c.cmd if @verbose
          end
          @data[i][:current] = outfile
          @data[j][:current] = outfile
          # run interleave cmd
        end
      end
      pair = " -p "
    else
      raise RuntimeError.new("single reads not currently supported as input")
    end
    # unpaired reads - cat all the single unpaired reads together
    cat_cmd = "cat "
    found=false
    @data.each_with_index do |a, i|
      if a[:unpaired]
        file = a[:unpaired]
        cat_cmd << " #{file} "
        found=true
      end
    end
    single_output = "#{@output_dir}/single_reads.fq"
    cat_cmd << " > #{single_output}"
    if !File.exist?(single_output) and found
      cat = Cmd.new(cat_cmd)
      cat.run
      puts cat.cmd if @verbose
    end
    # khmer
    set = Set.new
    @data.each_with_index do |a, i|
      set << a[:current]
    end
    file_list = set.to_a.join(" ")
    x = (@memory*1e9/buckets).to_i
    outfile = "#{@output_dir}/khmered.fq"

    cmd = "#{@khmer} #{pair} --quiet "
    cmd << " --ksize #{kmer} --cutoff #{cutoff} "
    cmd << " --n_tables #{buckets}  --min-tablesize #{x} "
    cmd << " --fault-tolerant --out #{outfile} "
    cmd << "#{file_list}"
    if !File.exist?(outfile)
      c = Cmd.new(cmd)
      c.run
      puts c.cmd if @verbose
    end

    if found
      # khmer on single unpaired reads
      outfile_single = "#{@output_dir}/#{@data[0][:type]}_khmered_single.fq"
      cmd_single = "#{@khmer} -q -k #{kmer} -C #{cutoff} -f -o #{outfile_single}"
      cmd_single << " #{single_output}"
      if !File.exist?(outfile_single)
        c = Cmd.new(cmd_single)
        c.run
        puts c.cmd if @verbose
      end
    end

    #deinterleave paired reads
    khmer_left = "#{@output_dir}/#{@data[0][:type]}.left.fq"
    khmer_right = "#{@output_dir}/#{@data[1][:type]}.right.fq"
    self.deinterleave(outfile, khmer_left, khmer_right)

    @data.each_with_index.each_slice(2) do |(a,i), (b,j)|
      @data[i][:current] = khmer_left
      @data[j][:current] = khmer_right
    end
  end

  def deinterleave(file, output_left, output_right)
    #puts "deinterleaving #{file} to make #{output_left} and #{output_right}"
    raise RuntimeError.new("Can't find #{file}") if !File.exist?(file)
    fastq = File.open(file)
    left = File.open("#{output_left}", "w")
    right= File.open("#{output_right}", "w")

    name1 = fastq.readline rescue nil
    seq1 = fastq.readline rescue nil
    plus1 = fastq.readline rescue nil
    qual1 = fastq.readline rescue nil
    name2 = fastq.readline rescue nil
    seq2 = fastq.readline rescue nil
    plus2 = fastq.readline rescue nil
    qual2 = fastq.readline rescue nil

    while name1 != nil and name2 != nil
      left.write(name1)
      left.write(seq1)
      left.write(plus1)
      left.write(qual1)
      right.write(name2)
      right.write(seq2)
      right.write(plus2)
      right.write(qual2)

      name1 = fastq.readline rescue nil
      seq1 = fastq.readline rescue nil
      plus1 = fastq.readline rescue nil
      qual1 = fastq.readline rescue nil
      name2 = fastq.readline rescue nil
      seq2 = fastq.readline rescue nil
      plus2 = fastq.readline rescue nil
      qual2 = fastq.readline rescue nil
    end
    fastq.close
    left.close
    right.close
  end

  def bbnorm(k=31, target_coverage=20, bits=8, tables=3)

    gem_dir = Gem.loaded_specs['preprocessor'].full_gem_path
    bbnorm = File.join(gem_dir, 'bin', 'bbmap', 'bbnorm.sh')
    # download
    if !File.exist?(bbnorm)
      dl = "#{gem_dir}/bin/bbmap.tar.gz"
      cmd = "wget http://kent.dl.sourceforge.net/project/bbmap/"
      cmd << "BBMap_33.08_java7.tar.gz -O #{dl}"

      wget = Cmd.new(cmd)
      wget.run
      puts wget.stdout if @verbose
      puts wget.cmd if @verbose
      if wget.status.success?
        cmd = "tar xzf #{dl} --directory "
        cmd << "#{gem_dir}/bin"
        untar = Cmd.new(cmd)
        untar.run
        puts untar.stdout if @verbose
        puts untar.cmd if @verbose
        if !untar.status.success?
          raise RuntimeError.new("untar failed")
        end
      else
        raise RuntimeError.new("download failed")
      end
    end
    # run for each pair
    @data.each_with_index.each_slice(2) do |(a,i), (b,j)|
      left  = a[:current]
      right = b[:current]
      single_left = a[:unpaired] if a[:unpaired]
      single_right = b[:unpaired] if b[:unpaired]
      outfile = "#{left}.keep" # interleaved output of left and right
      if left and right
        # parameters
        lowthresh=1 # Kmers at this and below are always considered errors
        mindepth=1  # Kmers with depth below this number will not be included
                    #   when calculating the depth of a read
        minkmers=15 # Reads must have at least this many kmers over min depth
                    #   to be retained.
        cmd = "#{bbnorm} in=#{left} in2=#{right} "
        cmd << "out=#{outfile} "
        # cmd << "outt=#{toss} " # not sure if this is necessary
        cmd << "k=#{k} "
        cmd << "hashes=#{tables} "
        cmd << "bits=#{bits} "
        cmd << "lowthresh=#{lowthresh} "
        cmd << "mindepth=#{mindepth} "
        cmd << "minkmers=#{mindepth} "
        cmd << "target=#{target_coverage} "
        cmd << "threads=#{@threads}"
        norm = Cmd.new(cmd)
        puts norm.cmd if @verbose
        norm.run
        if !norm.status.success?
          puts norm.stdout
          puts norm.stderr
          raise RuntimeError.new("something went wrong with bbnorm")
        end
        # now deinterleave output files back into left and right

        leftoutput = "#{@output_dir}/"
        leftoutput << "#{a[:type]}-#{a[:rep]}-#{a[:pair]}.bbnorm.fq"
        rightoutput = "#{@output_dir}/"
        rightoutput << "#{b[:type]}-#{b[:rep]}-#{b[:pair]}.bbnorm.fq"
        deinterleave(outfile, leftoutput, rightoutput)
        # delete intermediate
        File.delete(outfile)
        # save names
        a[:current] = leftoutput
        b[:current] = rightoutput
      end

      if single_left
        single_left_output = "#{@output_dir}/"
        single_left_output << "#{a[:type]}-#{a[:rep]}-#{a[:pair]}U.bbnorm.fq"
        cmd = "#{bbnorm} in=#{single_left} "
        cmd << "out=#{single_left_output} "
        cmd << "k=#{k} "
        cmd << "hashes=#{tables} "
        cmd << "bits=#{bits} "
        cmd << "lowthresh=#{lowthresh} "
        cmd << "mindepth=#{mindepth} "
        cmd << "target=#{target_coverage} "
        cmd << "threads=#{@threads}"
        norm = Cmd.new(cmd)
        puts norm.cmd if @verbose
        norm.run
        if !norm.status.success?
          puts norm.stdout
          puts norm.stderr
          raise RuntimeError.new("something went wrong with bbnorm")
        end
        @data[i][:unpaired] = single_left_output
      end

      if single_right
        single_right_output = "#{@output_dir}/"
        single_right_output << "#{b[:type]}-#{b[:rep]}-#{b[:pair]}U.bbnorm.fq"
        cmd = "#{bbnorm} in=#{single_right} "
        cmd << "out=#{single_right_output} "
        cmd << "k=#{k} "
        cmd << "hashes=#{tables} "
        cmd << "bits=#{bits} "
        cmd << "lowthresh=#{lowthresh} "
        cmd << "mindepth=#{mindepth} "
        cmd << "target=#{target_coverage} "
        cmd << "threads=#{@threads}"
        norm = Cmd.new(cmd)
        puts norm.cmd if @verbose
        norm.run
        if !norm.status.success?
          puts norm.stdout
          puts norm.stderr
          raise RuntimeError.new("something went wrong with bbnorm")
        end
        @data[j][:unpaired] = single_right_output
      end
    end
  end

  def norm

  end

  def get_output
    left_set = Set.new
    right_set = Set.new
    @data.each_with_index.each_slice(2) do |(a,i), (b,j)|
      left_set << @data[i][:current]
      right_set << @data[j][:current]
    end
    [left_set.to_a, right_set.to_a]
  end

  def run
    trim
    hammer
    khmer
  end
end
