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

  def initialize(input, output, verbose)
    @input = input
    @verbose = verbose
    @output_dir = File.expand_path(output)
    @trim_jar = "bin/trimmomatic-0.32.jar"
    @khmer = "normalize-by-median.py"
    @hammer_path = "/home/chris/documents/apps/SPAdes-3.0.0-Linux/bin/spades.py"
    @memory = 4 # TODO set this value
    @threads = 8 # TODO set this value
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
    phred = -1
    if max == 74 or max == 73
      phred = 33
    elsif max == 104 or max == 103
      phred = 64
    end
    return phred
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
        outfile_right = "#{@output_dir}/#{a[:type]}_#{a[:rep]}-#{b[:pair]}.t.fq"
        outfileU_left = "#{@output_dir}/#{a[:type]}_#{a[:rep]}-#{a[:pair]}.tU.fq"
        outfileU_right = "#{@output_dir}/#{a[:type]}_#{a[:rep]}-#{b[:pair]}.tU.fq"
        trim_cmd = "java -jar #{@trim_jar} PE "
        trim_cmd << " -phred#{self.detect_phred} "
        trim_cmd << " -threads 1 "
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
          puts trim_cmd if @verbose
          # run trim_cmd
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
    puts yaml
    File.open("#{@output_dir}/dataset.yaml","w") {|io| io.write yaml}
    cmd = "python #{@hammer_path} --dataset #{@output_dir}/dataset.yaml "
    cmd << " --only-error-correction "
    cmd << " --disable-gzip-output -m #{@memory} -t #{@threads}"
    cmd << " -o #{@output_dir}/hammer"
    puts cmd
    # loads output file location into @data[:current]
  end

  def khmer(kmer=23, cutoff=20, buckets=4, memory=8)
    x = (memory/buckets*1e9).to_i
    # interleave the input files if paired
    if @paired==2
      @data.each_with_index.each_slice(2) do |(a,i), (b,j)|
        input_left = a[:current]
        input_right = b[:current]
        if input_left and input_right
          outfile = "#{@output_dir}/#{File.basename(a[:current],"fq")}in.fq"
          cmd = "paste #{input_left} #{input_right} | paste - - - - | "
          cmd << "awk -v FS=\"\\t\" -v OFS=\"\\n\" \'{print(\"@read\"NR\":1\",$3,"
          cmd << "$5,$7,\"@read\"NR\":2\",$4,$6,$8)}\' > "
          cmd << " #{outfile}"
          puts cmd
          @data[i][:current] = outfile
          @data[j][:current] = outfile
          # run interleave cmd
        end
      end
      pair = " -p "
    end
    # khmer
    set = Set.new
    @data.each_with_index do |a, i|
      set << a[:current]
    end
    file_list = set.to_a.join(" ")
    outfile = "#{@output_dir}/khmered.fq"
    cmd = "#{@khmer} #{pair} -k #{kmer} -C #{cutoff} -f -o #{outfile} "
    cmd << "#{file_list}"
    puts cmd
  end
end
