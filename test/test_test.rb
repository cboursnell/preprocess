#!/usr/bin/env	ruby

require 'helper'
require 'tmpdir'

class TestPreprocessor < Test::Unit::TestCase

  context 'preprocessor' do

    setup do
      input = File.join(File.dirname(__FILE__), 'data', 'raw_data')
      # @output = Dir.mktmpdir
      @output = "dry_run"
      verbose = true
      threads = 1
      memory = 1
      @pre = Preprocessor.new(input, @output, verbose, threads, memory)
    end

    teardown do

    end

    should 'setup should run ok' do
      assert @pre
    end

    should 'trim reads using trimmomatic' do
      @pre.trim
      assert File.exist?("#{@output}/A_1-1.t.fq")
      assert File.exist?("#{@output}/A_1-2.t.fq")
      @pre.data.each do |hash|
        assert hash[:current]
        assert hash[:unpaired]
      end
    end

    should 'normalise reads with khmer' do
      @pre.khmer

      assert File.exist?("#{@output}/A.left.fq")
      assert File.exist?("#{@output}/A.right.fq")
      p @pre.data
    end

    should 'construct hammer input' do
      @pre.construct_hammer_input
      assert File.exist?("#{@output}/dataset.yaml")
    end

    should 'hammer reads' do
      @pre.hammer
    end

    should 'run bbnorm' do
      a = @pre.bbnorm
    end
  end
end
