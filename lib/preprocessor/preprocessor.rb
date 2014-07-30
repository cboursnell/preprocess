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
require 'which'
require 'bindeps'
include Which

module Preprocessor

  class MalformedInputError < StandardError
  end

  class Preprocessor

    attr_accessor :input, :output, :phred
    attr_reader :data

    def initialize(output, verbose, threads=1, memory=4)
      @verbose = verbose
      @filter = "http://zenodo.org/record/11091/files/rRNAplants.fa"
      @khmer = which("normalize-by-median.py").first
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

    def trimmomatic(minlen=40, windowsize=4, quality=15,
                    trailing=15, leading=15, mismatches=2)
      trimmer = Trimmomatic.new(@output_dir, @threads, minlen, windowsize,
                                quality, trailing, leading, mismatches)

      @data.each_with_index.each_slice(2) do |(left,i), (right,j)|
        trimmer.run(left, right)
      end

    end

    def skewer(end_quality=25, mean_quality=0, min_length=40)
      trimmer = Skewer.new(@output_dir, @theads, end_quality,
                           mean_quality, min_length)
      @data.each_with_index.each_slice(2) do |(left,i), (right,j)|
        trimmer.run(left, right)
      end

    end

    def hammer
      correcter = Hammer.new(@output_dir, @threads, @memory)
      @data.each_with_index.each_slice(2) do |(left,i), (right,j)|
        correcter.run(left, right)
      end
    end

    def khmer(kmer=23, cutoff=20, tables=4)
      # check that khmer is installed
      normalizer = Khmer.new(@output_dir, @threads, @memory, kmer, cutoff,
                            tables)
      @data.each_with_index.each_slice(2) do |(left,i), (right,j)|
        normalizer.interleave(left, right)
      end
      @data = normalizer.normalize
    end

    def bbnorm(k=31, target_coverage=20, bits=8, tables=3,
               lowthresh=1, mindepth=1, minkmers=15)
      normalizer = BBnorm.new(@output_dir, @threads, @memory,
              k, target_coverage, bits, tables, lowthresh, mindepth, minkmers)
      @data.each_with_index.each_slice(2) do |(left,i), (right,j)|
        normalizer.run(left, right)
      end
    end

    def facs(filter=nil, k=nil, false_positive=0.005, threshold=0.4)
      filter = @filter unless filter
      filterer = Facs.new(@output_dir, @threads, @memory, filter,
                          k, false_positive, threshold)

      @data.each_with_index.each_slice(2) do |(left,i), (right,j)|
        filterer.run(left, right)
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

end