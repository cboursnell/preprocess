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
  end
end
