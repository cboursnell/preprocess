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

class Preprocessor

  def initialize(input, verbose)
    @input = input
    @verbose = verbose
  end
end
