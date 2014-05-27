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
require 'trollop'

opts = Trollop::options do
  version "v0.1"
  banner <<-EOS
  mrna pipeline

  author: Chris Boursnell (cmb211@cam.ac.uk)
  EOS
  opt :input, "Input file describing fastq files",
      :required => true,
      :type => String
  opt :verbose, "Be verbose"
end

Trollop::die :input, "must exist" if !File.exist?(opts[:input]) if opts[:input]

class Preprocessor

  attr_accessor

  def initialize(input, verbose)
    @input = input
    @verbose = verbose
  end
end

p = Preprocessor.new(opts.input, opts.verbose)

