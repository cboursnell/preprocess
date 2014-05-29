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
require 'preprocessor'
require 'set'

class MalformedInputError < StandardError
end

class Preprocessor

  attr_accessor :input, :output, :phred

  def initialize(input, output, verbose, threads=1, memory=4)
    @input = input
    @verbose = verbose
    @output_dir = File.expand_path(output)
    @trim_jar = "bin/trimmomatic-0.32.jar"
    @khmer = "normalize-by-median.py"
    @hammer_path = "bin/SPAdes-3.0.0-Linux/bin/spades.py"
    @memory = memory
    @threads = threads
    @data = []
    if File.exist?(input)
      File.open("#{input}").each_line do |line|
        cols = line.chomp.split(",")
        if cols.size != 4
          raise MalformedInputError.new("Input file does not contain 4 columns")
        end
        raise RuntimeError.new("#{cols[0]} not found") if !File.exist?(cols[0])
        if cols[3].to_i != 1 and cols[3].to_i != 2
          raise RuntimeError.new("Pair should be 1 or 2") 
        end
        @data << { :file => File.expand_path(cols[0]),
                   :rep => cols[1].to_i,
                   :type => cols[2],
                   :pair => cols[3].to_i,
                   :current => File.expand_path(cols[0]) }
      end
    else
      raise RuntimeError, "#{input} does not exist"
    end
    @paired = @data.reduce(0) {|max,v| max=[max,v[:pair]].max}
  end

  def detect_phred
    file = File.open(@data[0][:file])
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
    if max == 74 or max == 73
      @phred = 33
    elsif max == 104 or max == 103
      @phred = 64
    end
    return @phred
  end

  def trim(minlen=40, windowsize=4, quality=15, trailing=15,
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
        trim_cmd << " #{a[:file]} #{b[:file]} "
        trim_cmd << " #{outfile_left} #{outfileU_left} "
        trim_cmd << "#{outfile_right} #{outfileU_right} "
        trim_cmd << " LEADING:#{leading} TRAILING:#{trailing} "
        trim_cmd << " SLIDINGWINDOW:#{windowsize}:#{quality} MINLEN:#{minlen}"
        @data[i][:current] = outfile_left
        @data[j][:current] = outfile_right
        @data[i][:unpaired] = outfileU_left
        @data[j][:unpaired] = outfileU_right
        if !File.exist?("#{outfile_left}")
          # puts trim_cmd if @verbose
          `#{trim_cmd}`
        else
          puts "trimmomatic already run on #{a[:file]}" if @verbose
        end
      end
    end
  end

  def hammer
    # wget http://spades.bioinf.spbau.ru/release3.0.0/SPAdes-3.0.0-Linux.tar.gz
    # tar -xzf SPAdes-3.0.0-Linux.tar.gz
    # cd SPAdes-3.0.0-Linux/bin/
    # --only-error-correction runs only read error correction (no assembly)
    # --disable-gzip-output forces error correction not to compress the 
    #                       corrected reads
    if !File.exist?(@hammer_path)
      wget_cmd = "wget http://spades.bioinf.spbau.ru/release3.0.0/"
      wget_cmd << "SPAdes-3.0.0-Linux.tar.gz -O bin/spades.tar.gz"
      `#{wget_cmd}`
      tar_cmd = "tar xzf bin/spades.tar.gz -C bin"
      `#{tar_cmd}`
    end
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
      yaml << "    ],\n  },\n  {\n"
      yaml << "    type: \"single\",\n    single reads: [\n"
      single.each_with_index do |single_read, j|
        yaml << "      \"#{single_read}\""
        yaml << "," if j < single.length-1
        yaml << "\n"
      end
      yaml << "    ]\n  }\n]\n"
    end
    # puts yaml
    File.open("#{@output_dir}/dataset.yaml","w") {|io| io.write yaml}
    cmd = "python #{@hammer_path} --dataset #{@output_dir}/dataset.yaml "
    cmd << " --only-error-correction "
    cmd << " --disable-gzip-output -m #{@memory} -t #{@threads}"
    cmd << " -o #{@output_dir}/hammer"
    puts cmd
    if !Dir.exist?("#{@output_dir}/hammer/corrected")
      `#{cmd}`
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
      leftU = a[:unpaired]
      rightU = b[:unpaired]
      left = File.basename(left).split(".")[0..-2].join(".")
      right = File.basename(right).split(".")[0..-2].join(".")
      leftU = File.basename(leftU).split(".")[0..-2].join(".")
      rightU = File.basename(rightU).split(".")[0..-2].join(".")
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
      corrected_single.each do |corr|
        if corr =~ /#{leftU}/
          @data[i][:unpaired] = corr
        end
        if corr =~ /#{rightU}/
          @data[j][:unpaired] = corr
        end
      end
    end
    # loads output file location into @data[:current]
  end

  def khmer(kmer=23, cutoff=20, buckets=4)
    x = (@memory/buckets*1e9).to_i
    # interleave the input files if paired
    pair = ""
    if @paired==2
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
          puts cmd
          if !File.exist?(outfile)
            `#{cmd}`
          end
          @data[i][:current] = outfile
          @data[j][:current] = outfile
          # run interleave cmd
        end
      end
      pair = " -p "
    end
    # unpaired reads
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
      `#{cat_cmd}`
    end
    # khmer
    set = Set.new
    @data.each_with_index do |a, i|
      set << a[:current]
    end
    file_list = set.to_a.join(" ")
    n = 4
    x = (@memory*1e9/n).to_i
    outfile = "#{@output_dir}/khmered.fq"
    cmd = "#{@khmer} #{pair} -q -k #{kmer} -C #{cutoff} -N #{n} -x #{x}"
    cmd << "-f -o #{outfile} "
    cmd << "#{file_list}"
    if !File.exist?(outfile)
      puts cmd
      `#{cmd}`
    end

    cmd_single = "#{@khmer} -q -k #{kmer} -C #{cutoff} -f -o #{outfile_single}"
    cmd_single << " #{single_output}"
    if !File.exist?(outfile_single)
      puts cmd_single
      `#{cmd_single}`
    end

    #deinterleave paired reads
    khmer_left = "#{@output_dir}/#{@data[0][:type]}.left.fq"
    khmer_right = "#{@output_dir}/#{@data[1][:type]}.right.fq"
    self.deinterleave(outfile, khmer_left , khmer_right)

    @data.each_with_index.each_slice(2) do |(a,i), (b,j)|
      @data[i][:current] = khmer_left
      @data[j][:current] = khmer_right
    end
  end

  def deinterleave(file, output_left, output_right)
    puts "deinterleaving #{file} to make #{output_left} and #{output_right}"
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
