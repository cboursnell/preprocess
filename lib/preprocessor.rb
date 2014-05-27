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

  def trim(minlen=40, windowsize=4, quality=15, trailing=15,
            leading=15, mismatches=2)
    if @paired==1
      @data.each_with_index do |a, i|
        puts "a = #{a}"
        puts "i = #{i}"
      end
    elsif @paired==2
      @data.each_with_index.each_slice(2) do |(a,i), (b,j)|
        puts "a = #{a}"
        puts "b = #{b}"
        puts "i = #{i}"
        puts "j = #{j}"
        outfile_left = "#{a[:cell]}_#{a[:rep]}-#{a[:pair]}.t.fq"
        outfile_right = "#{a[:cell]}_#{a[:rep]}-#{b[:pair]}.t.fq"
        outfileU_left = "#{a[:cell]}_#{a[:rep]}-#{a[:pair]}.tU.fq"
        outfileU_right = "#{a[:cell]}_#{a[:rep]}-#{b[:pair]}.tU.fq"
        trim_cmd = "java -jar #{@trim_jar} PE "
        trim_cmd << " -phred#{self.detect_phred} "
        trim_cmd << " -threads 1 "
        trim_cmd << " #{a[:file]} #{b[:file]} "
        trim_cmd << " #{outfile_left} #{outfileU_left} "
        trim_cmd << "#{outfile_right} #{outfileU_right} "
        trim_cmd << " LEADING:#{leading} TRAILING:#{trailing} "
        trim_cmd << " SLIDINGWINDOW:#{windowsize}:#{quality} MINLEN:#{minlen}"
        @data[i][:trimmed] = outfile_left
        @data[j][:trimmed] = outfile_right
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
end
