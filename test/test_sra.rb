#!/usr/bin/env	ruby

require 'helper'
require 'tmpdir'

class TestPreprocessor < Test::Unit::TestCase

  context 'preprocessor' do

    setup do
      @output = Dir.mktmpdir
      verbose = true
      threads = 1
      memory = 1
      @pre = Preprocessor::Preprocessor.new(@output, verbose, threads, memory)
    end

    teardown do
    end

    should 'decompress sra files into pairs' do
      @pre.data = [{:name => "sra_test",
                   :file => "test/data/SRR453569.sra",
                   :rep => 1,
                   :pair => 1,
                   :type => "yeast",
                   :current => "test/data/SRR453569.sra",
                   :processed => {} }]
      @pre.fastq_dump
      data = @pre.data
      assert File.exist?("#{@output}/SRR453569_1.fastq")
      assert File.exist?("#{@output}/SRR453569_2.fastq")
      assert_equal 2, data.size, "data size"
      assert_equal 1, data[0][:pair]
      assert_equal 2, data[1][:pair]
    end

  end
end
