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

class Preprocessor

  attr_accessor

  def initialize(input, verbose)
    @input = input
    @verbose = verbose
  end
end
