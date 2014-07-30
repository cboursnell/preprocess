require 'bindeps'

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
      gem_dir = Gem.loaded_specs['preprocessor'].full_gem_path
      @path = File.join(gem_dir, 'bin', 'trimmomatic-0.32.jar')
    end

    def run left, right=nil
      if right # paired
        outfile_left = "#{@outdir}/#{left[:type]}_#{left[:rep]}-#{left[:pair]}.t.fq"
        outfileU_left = "#{@outdir}/#{left[:type]}_#{left[:rep]}-#{left[:pair]}.tU.fq"

        outfile_right = "#{@outdir}/#{right[:type]}_#{right[:rep]}-#{right[:pair]}.t.fq"
        outfileU_right = "#{@outdir}/#{right[:type]}_#{right[:rep]}-#{right[:pair]}.tU.fq"

        trim_cmd = "java -jar #{@path} PE "
        trim_cmd << " -phred#{self.detect_phred left[:current]}"
        trim_cmd << " -threads #{@threads}"
        trim_cmd << " #{left[:current]} #{right[:current]}"
        trim_cmd << " #{outfile_left} #{outfileU_left}"
        trim_cmd << " #{outfile_right} #{outfileU_right}"
        trim_cmd << " LEADING:#{@leading} TRAILING:#{@trailing}"
        trim_cmd << " SLIDINGWINDOW:#{@windowsize}:#{@quality}"
        trim_cmd << " MINLEN:#{@minlen}"

        left[:current] = outfile_left
        left[:unpaired] = outfileU_left
        right[:current] = outfile_right
        right[:unpaired] = outfileU_right

        if !File.exist?("#{outfile_left}")
          cmd = Cmd.new(trim_cmd)
          cmd.run
        else
          puts "trimmomatic already run on #{left[:file]}"
        end

      else # unpaired

      end
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