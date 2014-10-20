require 'bindeps'
require 'which'
include Which

module Preprocessor

  class Trimmomatic

    def initialize(outdir, threads, minlen, windowsize, quality, trailing,
                   leading, mismatches)
      @outdir = outdir
      @threads = threads
      @minlen = minlen
      @windowsize = windowsize
      @quality = quality
      @trailing = trailing
      @leading = leading
      @mismatches = mismatches

      gem_dir = Gem.loaded_specs['preprocessor'].full_gem_path
      gem_deps = File.join(gem_dir, 'deps', 'trimmomatic.yaml')
      Bindeps.require gem_deps
      paths = which('trimmomatic-0.32.jar')
      if paths.empty?
        @path = File.join(ENV['GEM_HOME'], 'bin', 'trimmomatic-0.32.jar')
      else
        @path = paths.first
      end
      @java = which('java').first
      raise RuntimeError.new("java not installed") if !@java
    end

    def run left, right=nil
      cmd =""
      if right # paired
        outfile_left = "#{@outdir}/#{left[:type]}_#{left[:rep]}-#{left[:pair]}.t.fq"
        outfileU_left = "#{@outdir}/#{left[:type]}_#{left[:rep]}-#{left[:pair]}.tU.fq"

        outfile_right = "#{@outdir}/#{right[:type]}_#{right[:rep]}-#{right[:pair]}.t.fq"
        outfileU_right = "#{@outdir}/#{right[:type]}_#{right[:rep]}-#{right[:pair]}.tU.fq"

        cmd << "#{@java} -jar #{@path} PE "
        cmd << " -phred#{self.detect_phred left[:current]}"
        cmd << " -threads #{@threads}"
        cmd << " #{left[:current]} #{right[:current]}"
        cmd << " #{outfile_left} #{outfileU_left}"
        cmd << " #{outfile_right} #{outfileU_right}"
        cmd << " LEADING:#{@leading} TRAILING:#{@trailing}"
        cmd << " SLIDINGWINDOW:#{@windowsize}:#{@quality}"
        cmd << " MINLEN:#{@minlen}"

        left[:current] = outfile_left
        left[:unpaired] = outfileU_left
        right[:current] = outfile_right
        right[:unpaired] = outfileU_right

      else # unpaired
        outfile_left = "#{@outdir}/#{left[:type]}_#{left[:rep]}.t.fq"
        cmd << "#{@java} -jar #{@path} SE "
        cmd << " -phred#{self.detect_phred left[:current]}"
        cmd << " -threads #{@threads}"
        cmd << " #{left[:current]}"
        cmd << " #{outfile_left}"
        cmd << " LEADING:#{@leading} TRAILING:#{@trailing}"
        cmd << " SLIDINGWINDOW:#{@windowsize}:#{@quality}"
        cmd << " MINLEN:#{@minlen}"
        left[:current] = outfile_left
      end
      trim_cmd = Cmd.new(cmd)
      trim_cmd.run
      if !trim_cmd.status.success?
        msg = "trimmomatic failed\n#{trim_cmd.stdout}\n#{trim_cmd.stderr}"
        raise RuntimeError.new(msg)
      end
    end

    def stats file
      @file_histo = Array.new(101,0) if !@file_histo
      file = File.open(file[:current])
      name = file.readline
      seq = file.readline
      plus = file.readline
      qual = file.readline
      while name
        @file_histo[seq.chomp.length] ||= 0
        @file_histo[seq.chomp.length] += 1
        name = file.readline rescue nil
        seq = file.readline rescue nil
        plus = file.readline rescue nil
        qual = file.readline rescue nil
      end
      file.close
    end

    def get_stats
      tot = @file_histo.reduce(0) { |sum, c| sum += c }
      str = ""
      @file_histo.each_with_index do |count, length|
        if count > 0
          str << "#{length}\t#{count}\t#{(100 * count / tot.to_f).round(2)}%\n"
        end
      end
      str
    end

    def detect_phred filename
      file = File.open(filename)
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
      if max <= 74
        phred = 33
      elsif max <= 105
        phred = 64
      else
        raise RuntimeError.new("Couldn't determine phred formatting")
      end
      file.close
      return phred
    end

  end

end
