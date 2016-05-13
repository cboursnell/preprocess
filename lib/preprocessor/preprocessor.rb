#!/usr/bin/ruby

## # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
## mrna pipeline
##
## can run any of these
##   trimmomatic
##   skewer
##   bayeshammer
##   khmer
##   bbnorm
##   bowtie2
##   bwa
##   snap
##   express
##   salmon
##   more to come!
##
## created: 2014-05-27 Chris Boursnell (cmb211@cam.ac.uk)
## last updated: 2015-02-03
##
## # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

require 'rubygems'
require 'set'
require 'json'
require 'zlib'
require 'fixwhich'
require 'bindeps'
require 'fileutils'

module Preprocessor

  class MalformedInputError < StandardError
  end

  class Preprocessor

    attr_accessor :input, :output, :phred, :data

    def initialize(output, verbose, threads=1, memory=4)
      @verbose = verbose
      @filter = "http://zenodo.org/record/11091/files/rRNAplants.fa"
      @khmer = Which.which("normalize-by-median.py").first
      @output_dir = output ? File.expand_path(output) : Dir.pwd
      FileUtils.mkdir_p(@output_dir) unless Dir.exist?(@output_dir)
      memory.nil? ? @memory = 8 : @memory = memory
      threads.nil? ? @threads = 1 : @threads = threads
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

    def trimmomatic(minlen=40, windowsize=4, quality=5,
                    trailing=15, leading=15, mismatches=2, avgqual=15)
      trimmer = Trimmomatic.new(@output_dir, @threads, minlen, windowsize,
                                quality, trailing, leading, mismatches, avgqual)

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
      # don't get stats... takes too long
      # @data.each do |file|
      #   trimmer.stats file
      # end
      # File.open(File.join(@output_dir, "trimmomatic.stats"), "wb") do |out|
      #   out.write trimmer.get_stats
      # end
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
      stats = File.join(@output_dir, "hammer.stats")
      unless File.exist?(stats)
        @data.each do |info|
          correcter.stats info
        end
        File.open(stats, "wb") do |out|
          out.write correcter.get_stats
        end
      end
    end

    def hammer_batch
      correcter = Hammer.new(@output_dir, @threads, @memory)
      left_files = []
      right_files = []
      @data.each_with_index.each_slice(2) do |(left,i), (right,j)|
        left_files << left
        right_files << right
      end
      correcter.run_batch left_files, right_files

      stats = File.join(@output_dir, "hammer.stats")
      unless File.exist?(stats)
        @data.each do |info|
          correcter.stats info
        end
        File.open(stats, "wb") do |out|
          out.write correcter.get_stats
        end
      end
    end

    def rcorrector
      correcter = Rcorrector.new(@output_dir, @threads, @memory)
      if @paired==1
        @data.each_with_index do |left,i|
          left_files << left
        end
        correcter.run_single left_files
      elsif @paired==2
        left_files = []
        right_files = []
        @data.each_with_index.each_slice(2) do |(left,i), (right,j)|
          left[:processed][:correction] = "rcorrector"
          right[:processed][:correction] = "rcorrector"
          left_files << left
          right_files << right
        end
        print "Correcting with Rcorrector..." if @verbose
        correcter.run left_files, right_files
        puts " Done" if @verbose
        # TODO set output files into @data
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

    def batch
      cmd1 = "cat "
      cmd2 = "cat "
      @data.each_with_index.each_slice(2) do |(left,i), (right,j)|
        unless File.exist?(left[:current])
          abort "Can't find #{left[:current]}"
        end
        unless File.exist?(right[:current])
          abort "Can't find #{right[:current]}"
        end
        cmd1 << "#{left[:current]} "
        cmd2 << "#{right[:current]} "
      end
      left_output = "#{@output_dir}/#{@data[0][:name]}_1.fq"
      right_output = "#{@output_dir}/#{@data[1][:name]}_2.fq"
      cmd1 << " > #{left_output}"
      cmd2 << " > #{right_output}"
      cat1 = Cmd.new cmd1
      cat2 = Cmd.new cmd2
      unless File.exist?(left_output)
        puts "concatenating left files..." if @verbose
        cat1.run left_output
        unless cat1.status.success?
          puts "Cat1 failed"
        end
      end
      unless File.exist?(right_output)
        puts "concatenating right files..." if @verbose
        cat2.run right_output
        unless cat2.status.success?
          puts "Cat2 failed"
        end
      end
      @data[0][:current] = left_output
      @data[1][:current] = right_output
      @data = [@data[0], @data[1]]
    end

    def bbnorm(k=31, target_coverage=20, bits=16, tables=3,
               lowthresh=1, mindepth=1, minkmers=10)
      normaliser = BBnorm.new(@output_dir, @threads, @memory,
              k, target_coverage, bits, tables, lowthresh, mindepth, minkmers)
      @data.each_with_index.each_slice(2) do |(left,i), (right,j)|
        puts "bbnorm on #{left[:current]}..." if @verbose
        normaliser.run(left, right)
        left[:processed][:normalise] = "bbnorm"
        right[:processed][:normalise] = "bbnorm"
        File.open("#{@output_dir}/log", "wb") do |f|
          f.write(JSON.pretty_generate(@data))
        end
      end
    end

    def bowtie2(reference, expression=false)
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

    def filter(fasta)
      filter = Filter.new(@output_dir, @threads, fasta)
      if @paired==1
        @data.each_with_index do |left, i|
          puts "Filtering #{left[:type]}-#{left[:rep]}..." if @verbose
          filter.run(left)
          left[:processed][:align] = "filter"
        end
      elsif @paired==2
        @data.each_with_index.each_slice(2) do |(left,i), (right,j)|
          puts "Filtering #{left[:type]}-#{left[:rep]}..." if @verbose
          filter.run(left, right)
          left[:processed][:align] = "filter"
          right[:processed][:align] = "filter"
          File.open("#{@output_dir}/log", "wb") do |f|
            f.write(JSON.pretty_generate(@data))
          end
        end
      end
    end

    def bwa(reference, expression=false)
      aligner = Bwa.new(@output_dir, @threads, reference, expression)
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
        end
      end
      File.open("#{@output_dir}/log", "wb") do |f|
        f.write(JSON.pretty_generate(@data))
      end
    end

    def snap(reference, expression=false)
      aligner = Snap.new(@output_dir, @threads, reference, expression)
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

    def salmon(reference)
      salmon = Salmon.new(@output_dir, @threads, reference)

      if @paired == 1
        @data.each_with_index do |left, i|
          puts "Quantifying #{left[:current]}" if @verbose
          salmon.run(left)
          left[:processed][:expression] = "salmon"
          File.open("#{@output_dir}/log", "wb") do |f|
            f.write(JSON.pretty_generate(@data))
          end
        end
        results = "#{@output_dir}/salmon_results"
        File.open(results, "wb") do |out|
          @data.each_with_index do |left,i|
            name = left[:name]
            type = left[:type]
            rep = left[:rep]
            file = left[:salmon]
            out.write("#{name},#{type},#{rep},#{file}\n")
          end
        end
      else # paired == 2
        @data.each_with_index.each_slice(2) do |(left,i), (right,j)|
          puts "Quantifying #{left[:current]}" if @verbose
          salmon.run(left, right)
          left[:processed][:expression] = "salmon"
          right[:processed][:expression] = "salmon"
          File.open("#{@output_dir}/log", "wb") do |f|
            f.write(JSON.pretty_generate(@data))
          end
        end
        results = "#{@output_dir}/salmon_results"
        File.open(results, "wb") do |out|
          @data.each_with_index.each_slice(2) do |(left,i), (right,j)|
            name = left[:name]
            type = left[:type]
            rep = left[:rep]
            file = left[:salmon]
            out.write("#{name},#{type},#{rep},#{file}\n")
          end
        end
      end
    end

    # :name=>"test", :rep=>1, :type=>"A", :pair=>1,
    def combine_salmon
      @tpm = {}
      @counts = {}
      if @paired == 1
        @data.each_with_index do |left,i|
          name = "#{left[:name]}-#{left[:type]}-#{left[:rep]}"
          File.open(left[:salmon]).each do |line|
            unless line.start_with?("#")
              cols = line.chomp.split("\t")
              transcript = cols[0]
              tpm = cols[2].to_f
              counts = cols[3].to_f
              @tpm[transcript] ||= {}
              @tpm[transcript][name] = tpm
              @counts[transcript] ||= {}
              @counts[transcript][name] = counts
            end
          end
        end
      elsif @paired == 2
        @data.each_with_index.each_slice(2) do |(left,i), (right,j)|
          name = "#{left[:name]}-#{left[:type]}-#{left[:rep]}"
          File.open(left[:salmon]).each do |line|
            unless line.start_with?("#")
              cols = line.chomp.split("\t")
              transcript = cols[0]
              tpm = cols[2].to_f
              counts = cols[3].to_f
              @tpm[transcript] ||= {}
              @tpm[transcript][name] = tpm
              @counts[transcript] ||= {}
              @counts[transcript][name] = counts
            end
          end
        end
      end
      File.open("#{@output_dir}/combined_tpm.csv", "wb") do |out|
        out.write "transcript"
        @tpm[@tpm.keys.first].each do |file, tpm|
          out.write "\t#{file}"
        end
        out.write "\n"
        @tpm.each do |transcript, hash|
          out.write "#{transcript}"
          hash.each do |file, tpm|
            out.write "\t#{tpm}"
          end
          out.write "\n"
        end
      end
      File.open("#{@output_dir}/combined_counts.csv", "wb") do |out|
        out.write "transcript"
        @counts[@counts.keys.first].each do |file, counts|
          out.write "\t#{file}"
        end
        out.write "\n"
        @counts.each do |transcript, hash|
          out.write "#{transcript}"
          hash.each do |file, counts|
            out.write "\t#{counts}"
          end
          out.write "\n"
        end
      end
    end # combine_salmon

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
