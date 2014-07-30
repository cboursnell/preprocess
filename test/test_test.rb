#!/usr/bin/env	ruby

require 'helper'
require 'tmpdir'

class TestPreprocessor < Test::Unit::TestCase

  context 'preprocessor' do

    setup do
      input = File.join(File.dirname(__FILE__), 'data', 'raw_data')
      @output = Dir.mktmpdir
      verbose = true
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
      cmd = "rm #{bindir}/facs"
      `#{cmd}` if File.exist?(File.join(bindir, "facs"))
      cmd = "rm -rf #{bindir}/bbmap"
      `#{cmd}` if Dir.exist?(File.join(bindir, "bbmap"))
      cmd = "rm #{bindir}/bbmap.tar.gz"
      `#{cmd}` if File.exist?(File.join(bindir, "bbmap.tar.gz"))
      cmd = "rm -rf #{bindir}/SPAdes-3.1.0-Linux"
      `#{cmd}` if Dir.exist?(File.join(bindir, "SPAdes-3.1.0-Linux"))
      cmd = "rm #{bindir}/spades.tar.gz"
      `#{cmd}` if File.exist?(File.join(bindir, "spades.tar.gz"))
    end

    should 'setup should run ok' do
      assert @pre
    end

    should 'install trimmomatic' do
      @pre.trimmomatic
      bindir = File.join(ENV['GEM_HOME'], 'bin')
      cmd = "java -jar #{bindir}/trimmomatic-0.32.jar"
      java_cmd = Preprocessor::Cmd.new(cmd)
      java_cmd.run
      assert_equal 357, java_cmd.stderr.length
      assert_equal "Usage: ", java_cmd.stderr.split("\n").first
    end

    should 'trim reads using trimmomatic' do
      @pre.trimmomatic
      assert File.exist?("#{@output}/A_1-1.t.fq"), "file doesn't exist"
      assert File.exist?("#{@output}/A_1-2.t.fq"), "file doesn't exist"
      @pre.data.each do |hash|
        assert hash[:current]
        assert hash[:unpaired]
      end
    end

    should 'install skewer' do
      @pre.skewer
      cmd = "skewer-0.1.117-linux-x86_64 --help"
      skewer_cmd = Preprocessor::Cmd.new(cmd)
      skewer_cmd.run
      stdout = skewer_cmd.stdout.split("\n")
      str="Skewer (A fast and accurate adapter trimmer for paired-end reads)"
      ver="Version 0.1.117 (updated in July 12, 2014), Author: Hongshan Jiang"
      assert_equal str, stdout[0]
      assert_equal ver, stdout[1]
    end

    should 'trim reads using skewer' do
      @pre.skewer
      @pre.data.each do |hash|
        assert File.exist?(hash[:current]), "file exist"
        assert hash[:current]=~/test-A-[12]-pair/, "file name"
      end
    end

    should 'load reads without input file' do
      verbose = false
      threads = 1
      memory = 1
      left = File.join(File.dirname(__FILE__), 'data', 'A-1-1.fq')
      left << ",#{File.join(File.dirname(__FILE__), 'data', 'A-2-1.fq')}"
      right = File.join(File.dirname(__FILE__), 'data', 'A-1-2.fq')
      right << ",#{File.join(File.dirname(__FILE__), 'data', 'A-2-2.fq')}"
      pre = Preprocessor::Preprocessor.new(@output, verbose, threads, memory)
      pre.load_reads(left, right, "A")
      pre.trimmomatic
      pre.hammer
      pre.bbnorm
    end

    should "gunzip files and leave original in place" do
      left = File.join(File.dirname(__FILE__), 'data', 'A-1-1.fq.gz')
      right = File.join(File.dirname(__FILE__), 'data', 'A-1-2.fq.gz')
      output = Dir.mktmpdir
      pre = Preprocessor::Preprocessor.new(output, false, 1, 1)
      pre.load_reads(left, right, "A")
      pre.gunzip
      gz1 = File.join(File.dirname(__FILE__), 'data', 'A-1-1.fq.gz')
      gz2 = File.join(File.dirname(__FILE__), 'data', 'A-1-2.fq.gz')
      assert File.exist?(gz1), "original file 1 exists"
      assert File.exist?(gz2), "original file 2 exists"
      assert File.exist?("#{output}/A-1-1.fq"), "gunzipped file 1 exists"
      assert File.exist?("#{output}/A-1-2.fq"), "gunzipped file 2 exists"
      assert_equal "#{output}/A-1-1.fq", pre.data[0][:current]
      assert_equal "#{output}/A-1-2.fq", pre.data[1][:current]
    end

    should 'normalise trimmed reads with khmer' do
      @pre.trimmomatic
      @pre.khmer(23, 5, 4)
      assert File.exist?("#{@output}/test.khmered.fq-left.fq")
      assert File.exist?("#{@output}/test.khmered.fq-right.fq")
      assert File.exist?("#{@output}/test.unpaired.fq")
    end

    should 'hammer reads' do
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
      @pre.trimmomatic
      @pre.hammer
      files = []
      files << "#{@output}/hammer-test-A-1/corrected/A_1-1.t.00.0_0.cor.fastq"
      files << "#{@output}/hammer-test-A-1/corrected/A_1-2.t.00.0_0.cor.fastq"
      files << "#{@output}/hammer-test-A-1/left_unpaired.fq"
      files << "#{@output}/hammer-test-A-1/right_unpaired.fq"
      files << "#{@output}/hammer-test-A-2/corrected/A_2-1.t.00.0_0.cor.fastq"
      files << "#{@output}/hammer-test-A-2/corrected/A_2-2.t.00.0_0.cor.fastq"
      files << "#{@output}/hammer-test-A-2/left_unpaired.fq"
      files << "#{@output}/hammer-test-A-2/right_unpaired.fq"
      files.each do |f|
        assert File.exist?(f), "file #{f} not found"
      end
    end

    should 'run bbnorm' do
      @pre.bbnorm
      assert File.exist?("#{@output}/test-A-1-1.bbnorm.fq")
      assert File.exist?("#{@output}/test-A-1-2.bbnorm.fq")
      assert File.exist?("#{@output}/test-A-2-1.bbnorm.fq")
      assert File.exist?("#{@output}/test-A-2-2.bbnorm.fq")
    end

    should 'trim and then bbnorm' do
      @pre.trimmomatic
      @pre.bbnorm
      assert File.exist?("#{@output}/test-A-1-1.bbnorm.fq")
      assert File.exist?("#{@output}/test-A-1-2.bbnorm.fq")
      assert File.exist?("#{@output}/test-A-2-1.bbnorm.fq")
      assert File.exist?("#{@output}/test-A-2-2.bbnorm.fq")
    end

    should 'filter reads with facs' do
      @pre.facs
      files = []
      files << "#{@output}/rRNAplants.fa.bloom"
      files << "#{@output}/A-1-1_rRNAplants.fa_clean.fastq"
      files << "#{@output}/A-1-2_rRNAplants.fa_clean.fastq"
      files << "#{@output}/A-2-1_rRNAplants.fa_clean.fastq"
      files << "#{@output}/A-2-2_rRNAplants.fa_clean.fastq"
      files.each do |f|
        assert File.exist?(f), "file #{f} not found"
      end
      assert @pre.data[0][:current] == files[1], "data is correct"
    end

    should 'filter reads with facs and specified file' do
      file = File.join(File.dirname(__FILE__), 'data', 'filter.fasta')
      @pre.facs(file)
      files = []
      files << "#{@output}/filter.fasta.bloom"
      files << "#{@output}/A-1-1_filter.fasta_clean.fastq"
      files << "#{@output}/A-1-2_filter.fasta_clean.fastq"
      files << "#{@output}/A-2-1_filter.fasta_clean.fastq"
      files << "#{@output}/A-2-2_filter.fasta_clean.fastq"
      files.each do |f|
        assert File.exist?(f), "file #{f} not found"
      end
      assert @pre.data[0][:current] == files[1], "data is correct"
    end

    should 'run my normaliser' do
      @pre.norm
    end
  end
end
