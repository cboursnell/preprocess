#!/usr/bin/env  ruby

require 'helper'
require 'tmpdir'

class TestQuantify < Test::Unit::TestCase

  context 'quantification' do

    setup do
      input = File.join(File.dirname(__FILE__), 'data', 'raw_data')
      @output = Dir.mktmpdir
      verbose = false
      threads = 1
      memory = 1
      @pre = Preprocessor::Preprocessor.new(@output, verbose, threads, memory)
      @pre.load_input(input)
    end

    teardown do
      # delete output folder
      cmd = "rm -rf #{@output}"
      `#{cmd}`
    end

    # should 'run express' do
    #   reference = File.join(File.dirname(__FILE__), 'data', 'reference.fa')
    #   @pre.bowtie2(reference, true)
    # end

    should 'run salmon' do
      puts "output tmp dir: #{@output}"
      reference = File.join(File.dirname(__FILE__), 'data', 'reference.fa')
      @pre.bowtie2(reference, true)
      @pre.salmon(reference)
      assert File.exist?("#{@output}/quant.sf")
    end

  end
end
