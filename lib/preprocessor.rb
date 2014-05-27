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

class MalformedInputError < StandardError
end

class Preprocessor

  def initialize(input, output, verbose)
    @input = input
    @verbose = verbose
    @trim_jar = "bin/trimmomatic-0.32.jar"
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
        @data << { :file => cols[0],
                   :rep => cols[1].to_i,
                   :type => cols[2],
                   :pair => cols[3].to_i }
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
  end
end
