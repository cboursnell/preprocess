#!/usr/bin/ruby

## # # # # # # # # # # # # #
## mrna pipeline
##
## can run any of these
##   trimmomatic
##   skewer
##   bayeshammer
##   khmer
##   bbnorm
##   more to come!
##
## created: 2014-05-27 Chris Boursnell (cmb211@cam.ac.uk)
##
## # # # # # # # # # # # # #

require 'rubygems'
require 'set'
require 'json'
require 'zlib'
require 'which'
require 'bindeps'
include Which

module Preprocessor

  class MalformedInputError < StandardError
  end

  class Preprocessor

    attr_accessor :input, :output, :phred, :data

    def initialize(output, verbose, threads=1, memory=4)
      @verbose = verbose
      @filter = "http://zenodo.org/record/11091/files/rRNAplants.fa"
      @khmer = which("normalize-by-median.py").first
      @output_dir = output ? File.expand_path(output) : Dir.pwd
      Dir.mkdir(@output_dir) unless Dir.exist?(@output_dir)
      @memory = memory
      @threads = threads
      @data = []
    end

    def load_input(input)
      if File.exist?(input)
        File.open("#{input}").each_line do |line|
          cols = line.chomp.split(",")
          if cols.size != 5
            msg = "Input file does not contain 5 columns\n"
            msg << "Please refer to documentation or use --morehelp option"
            raise MalformedInputError.new(msg)
          end
          if !File.exist?(cols[1])
            raise RuntimeError.new("#{cols[1]} not found")
          end
          if cols[4].to_i != 1 and cols[4].to_i != 2
            raise RuntimeError.new("Pair should be 1 or 2")
          end
          @data << { :name => cols[0],
                     :file => File.expand_path(cols[1]),
                     :rep => cols[2].to_i,
                     :type => cols[3],
                     :pair => cols[4].to_i,
                     :current => File.expand_path(cols[1]),
                     :processed => {} }
        end
        @paired = @data.reduce(0) { |max,v| max=[max,v[:pair]].max }
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
                   :current => File.expand_path(a),
                   :processed => {} }
        # right
        @data << { :name => name,
                   :file => File.expand_path(b),
                   :rep => rep,
                   :type => name,
                   :pair => 2,
                   :current => File.expand_path(b),
                   :processed => {} }
        rep += 1
      end
      @paired = 2
    end

    def load_single_reads(left, name)
      @data = []
      rep = 1
      left.split(",").each do |a|
        @data << { :name => name,
                   :file => File.expand_path(a),
                   :rep => rep,
                   :type => name,
                   :pair => 1,
                   :current => File.expand_path(a),
                   :processed => {} }
        rep += 1
      end
      @paired = 1
    end

    def gunzip
      @data.each do |info|
        if info[:current]=~/\.gz$/
          output_filename = File.basename(info[:current].split(".gz").first)
          output_filename = File.join(@output_dir, output_filename)
          File.open("#{output_filename}", "wb") do |out|
            Zlib::GzipReader.open(info[:current]) do |gz|
              out.write(gz.read)
            end
          end
          info[:current] = output_filename
          info[:processed][:unzip] = "gunzip"
          File.open("#{@output_dir}/log", "wb")  do |f|
            f.write(JSON.pretty_generate(@data))
          end
        end
      end
    end

    def fastq_dump
      tmp = @data.dup
      tmp.each do |info|
        if info[:current]=~/\.sra$/
          gem_dir = Gem.loaded_specs['preprocessor'].full_gem_path
          gem_deps = File.join(gem_dir, 'deps', 'fastq-dump.yaml')
          Bindeps.require gem_deps
          cmd = "fastq-dump.2.4.0 --origfmt --split-3 #{info[:current]} "
          cmd << "--outdir #{@output_dir}"
          extract = Cmd.new(cmd)
          extract.run
          basename = File.basename(info[:current], File.extname(info[:current]))
          left = info.dup
          right = info.dup
          left[:current] = "#{@output_dir}/#{basename}_1.fastq"
          right[:current] = "#{@output_dir}/#{basename}_2.fastq"
          left[:pair] = 1
          right[:pair] = 2
          left[:processed][:uncompress] = "fastq-dump"
          right[:processed][:uncompress] = "fastq-dump"
          @data << left
          @data << right
        end
      end
    end

    def set_output(output_dir)
      @output_dir = File.expand_path(output_dir)
    end

    def trimmomatic(minlen=40, windowsize=4, quality=15,
                    trailing=15, leading=15, mismatches=2)
      trimmer = Trimmomatic.new(@output_dir, @threads, minlen, windowsize,
                                quality, trailing, leading, mismatches)

      if @paired == 1
        @data.each_with_index do |left, i|
          puts "Trimming #{left[:type]}-#{left[:rep]}..." if @verbose
          right=nil
          trimmer.run(left, right)
          left[:processed][:trim] = "trimmomatic"
          File.open("#{@output_dir}/log", "wb") do |f|
            f.write(JSON.pretty_generate(@data))
          end
        end
      elsif @paired == 2
        @data.each_with_index.each_slice(2) do |(left,i), (right,j)|
          puts "Trimming #{left[:type]}-#{left[:rep]}..." if @verbose
          trimmer.run(left, right)
          left[:processed][:trim] = "trimmomatic"
          right[:processed][:trim] = "trimmomatic"
          File.open("#{@output_dir}/log", "wb") do |f|
            f.write(JSON.pretty_generate(@data))
          end
        end
      end
      @data.each do |file|
        trimmer.stats file
      end
      File.open(File.join(@output_dir, "trimmomatic.stats"), "wb") do |out|
        out.write trimmer.get_stats
      end
    end

    def skewer(end_quality=25, mean_quality=0, min_length=40)
      trimmer = Skewer.new(@output_dir, @theads, end_quality,
                           mean_quality, min_length)
      @data.each_with_index.each_slice(2) do |(left,i), (right,j)|
        trimmer.run(left, right)
        left[:processed][:trim] = "skewer"
        right[:processed][:trim] = "skewer"
        File.open("#{@output_dir}/log", "wb") do |f|
          f.write(JSON.pretty_generate(@data))
        end
      end
    end

    def hammer
      correcter = Hammer.new(@output_dir, @threads, @memory)
      @data.each_with_index.each_slice(2) do |(left,i), (right,j)|
        puts "Correcting #{left[:type]}-#{left[:rep]}..." if @verbose
        correcter.run(left, right)
        left[:processed][:correction] = "bayeshammer"
        right[:processed][:correction] = "bayeshammer"
        File.open("#{@output_dir}/log", "wb") do |f|
          f.write(JSON.pretty_generate(@data))
        end
      end
      @data.each do |file|
        correcter.stats file
      end
      File.open(File.join(@output_dir, "hammer.stats"), "wb") do |out|
        out.write correcter.get_stats
      end

    end

    def khmer(kmer=23, cutoff=20, tables=4)
      # check that khmer is installed
      normaliser = Khmer.new(@output_dir, @threads, @memory, kmer, cutoff,
                            tables)
      @data.each_with_index.each_slice(2) do |(left,i), (right,j)|
        normaliser.interleave(left, right)
      end
      @data = normaliser.normalise
      @data.each do |file|
        file[:processed][:normalise] = "khmer"
        File.open("#{@output_dir}/log", "wb") do |f|
          f.write(JSON.pretty_generate(@data))
        end
      end
    end

    def bbnorm(k=31, target_coverage=20, bits=8, tables=3,
               lowthresh=1, mindepth=1, minkmers=15)
      normaliser = BBnorm.new(@output_dir, @threads, @memory,
              k, target_coverage, bits, tables, lowthresh, mindepth, minkmers)
      @data.each_with_index.each_slice(2) do |(left,i), (right,j)|
        normaliser.run(left, right)
        left[:processed][:normalise] = "bbnorm"
        right[:processed][:normalise] = "bbnorm"
        File.open("#{@output_dir}/log", "wb") do |f|
          f.write(JSON.pretty_generate(@data))
        end
      end
    end

    def bowtie2(reference, expression)
      aligner = Bowtie2.new(@output_dir, @threads, reference, expression)
      if @paired==1
        @data.each_with_index do |left, i|
          puts "Aligning #{left[:type]}-#{left[:rep]}..." if @verbose
          aligner.run(left)
          left[:processed][:align] = "bowtie2"
        end
      elsif @paired==2
        @data.each_with_index.each_slice(2) do |(left,i), (right,j)|
          puts "Aligning #{left[:type]}-#{left[:rep]}..." if @verbose
          aligner.run(left, right)
          left[:processed][:align] = "bowtie2"
          right[:processed][:align] = "bowtie2"
          File.open("#{@output_dir}/log", "wb") do |f|
            f.write(JSON.pretty_generate(@data))
          end
        end
      end
      File.open(File.join(@output_dir, "bowtie2.stats"), "wb") do |out|
        out.write aligner.get_stats
      end
    end

    def snap(reference)
      aligner = Snap.new(@output_dir, @threads, reference)
      @data.each_with_index.each_slice(2) do |(left,i), (right,j)|
        puts "Aligning #{left[:type]}-#{left[:rep]}..." if @verbose
        aligner.run(left, right)
        left[:processed][:align] = "snap"
        right[:processed][:align] = "snap"
        File.open("#{@output_dir}/log", "wb") do |f|
          f.write(JSON.pretty_generate(@data))
        end
      end
      File.open(File.join(@output_dir, "snap.stats"), "wb") do |out|
        out.write aligner.get_stats
      end
    end

    def express(reference)
      expression = Express.new(@output_dir, @threads, reference)
      @data.each_with_index.each_slice(2) do |(left,i), (right,j)|
        expression.run(left, right)
        left[:processed][:expression] = "eXpress"
        right[:processed][:expression] = "eXpress"
        File.open("#{@output_dir}/log", "wb") do |f|
          f.write(JSON.pretty_generate(@data))
        end
      end
      # output a csv of name,type,rep,express_file_path
      # to be passed into Rscript
      results = "#{@output_dir}/express_results"
      File.open(results, "wb") do |out|
        @data.each_with_index.each_slice(2) do |(left,i), (right,j)|
          name = left[:name]
          type = left[:type]
          rep = left[:rep]
          file = left[:express]
          out.write("#{name},#{type},#{rep},#{file}\n")
        end
      end
      #Usage: ebseq.R [-[-help|h]] [-[-threads|t] <integer>]
      # [-[-files|f] <character>] [-[-output|o] <character>]
      ### rscript
      # gem_dir = Gem.loaded_specs['preprocessor'].full_gem_path
      # ebseq_path = File.join(gem_dir, "lib", "ebseq.R")
      # cmd = "Rscript #{ebseq_path} -t #{@threads} "
      # cmd << " -f #{results} -o #{@output_dir}"

      # ebseq = Cmd.new(cmd)
      # ebseq.run
      # puts ebseq.stdout
      # puts ebseq.stderr
    end

    def cbnorm # yet to be implemented

    end

    def get_output
      left_set = Set.new
      right_set = Set.new
      single_set = Set.new
      @data.each_with_index.each_slice(2) do |(a,i), (b,j)|
        left_set << a[:current]
        right_set << b[:current]
        single_set << a[:unpaired] if a[:unpaired]
        single_set << b[:unpaired] if b[:unpaired]
      end
      [left_set.to_a, right_set.to_a, single_set.to_a]
    end

  end

end