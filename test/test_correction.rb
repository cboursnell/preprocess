#!/usr/bin/env	ruby

require 'helper'
require 'tmpdir'

class TestCorrection < Test::Unit::TestCase

  context 'correction' do

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
      gem_dir = Gem.loaded_specs['preprocessor'].full_gem_path
      bindir = File.join(gem_dir, "bin")
      cmd = "rm #{bindir}/trimmomatic-0.32.jar"
      `#{cmd}` if File.exist?(File.join(bindir, "trimmomatic-0.32.jar"))
      cmd = "rm -rf #{bindir}/SPAdes-3.5.0-Linux"
      `#{cmd}` if Dir.exist?(File.join(bindir, "SPAdes-3.5.0-Linux"))
      cmd = "rm #{bindir}/spades.tar.gz"
      `#{cmd}` if File.exist?(File.join(bindir, "spades.tar.gz"))
    end

    should 'hammer reads' do
      puts @output
      @pre.hammer
      files = []
      files << "#{@output}/hammer-test-A-1/corrected/A-1-1.00.0_0.cor.fastq"
      files << "#{@output}/hammer-test-A-1/corrected/A-1-2.00.0_0.cor.fastq"
      files << "#{@output}/hammer-test-A-1/left_unpaired.fq"
      files << "#{@output}/hammer-test-A-2/corrected/A-2-1.00.0_0.cor.fastq"
      files << "#{@output}/hammer-test-A-2/corrected/A-2-2.00.0_0.cor.fastq"
      files << "#{@output}/hammer-test-A-2/left_unpaired.fq"
      files.each do |f|
        assert File.exist?(f), "file #{f} not found"
      end
    end

    should 'hammer trimmed reads' do
      puts @output
      @pre.trimmomatic
      @pre.hammer
      files = []
      files << "#{@output}/hammer-test-A-1/corrected/test_A_1-1.t.00.0_0.cor.fastq"
      files << "#{@output}/hammer-test-A-1/corrected/test_A_1-2.t.00.0_0.cor.fastq"
      files << "#{@output}/hammer-test-A-1/left_unpaired.fq"
      files << "#{@output}/hammer-test-A-1/right_unpaired.fq"
      files << "#{@output}/hammer-test-A-2/corrected/test_A_2-1.t.00.0_0.cor.fastq"
      files << "#{@output}/hammer-test-A-2/corrected/test_A_2-2.t.00.0_0.cor.fastq"
      files << "#{@output}/hammer-test-A-2/left_unpaired.fq"
      files << "#{@output}/hammer-test-A-2/right_unpaired.fq"
      files.each do |f|
        assert File.exist?(f), "file #{f} not found"
      end
    end

  end
end
