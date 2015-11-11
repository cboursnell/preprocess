require 'bindeps'
require 'csv'
require 'fixwhich'

module Preprocessor

  class Filter

    def initialize(outdir, threads, fasta)
      @outdir = outdir
      @threads = threads
      @fasta = File.expand_path(fasta)

      gem_dir = Gem.loaded_specs['preprocessor'].full_gem_path
      gem_deps = File.join(gem_dir, 'deps', 'bowtie2.yaml')
      Bindeps.require gem_deps
      @build = Which.which('bowtie2-build').first
      @bowtie = Which.which('bowtie2').first
    end

    def build_index
      # construct index if it doesn't exist
      index = File.join(@outdir, File.basename(@fasta))
      index_file = "#{index}.1.bt2"
      if !File.exist?(index_file)
        cmd = "#{@build}"
        cmd << " #{@fasta}"
        cmd << " #{index}"
        build = Cmd.new(cmd)
        build.run index_file
        if !build.status.success?
          msg = "Something went wrong with bowtie2-build\n"
          msg << "#{cmd}\n"
          msg << "#{build.stdout}"
          msg << "#{build.stderr}"
          raise RuntimeError.new(msg)
        end
      end
      index
    end

    def run left, right=nil
      index = build_index
      l = File.basename(left[:current])
      r = File.basename(right[:current]) if right
      unal = "#{left[:name]}_#{left[:type]}-#{left[:rep]}-#{left[:pair]}.f.fq"
      unal = File.join(@outdir, unal)

      cmd = "#{@bowtie}"
      cmd << " -x #{index}"
      if right
        cmd << " -1 #{left[:current]}"
        cmd << " -2 #{right[:current]}"
      else
        cmd << " -U #{left[:current]}"
      end
      cmd << " -p #{@threads}"
      cmd << " --very-sensitive "
      if right
        cmd << " --un-conc #{unal}"
      else
        cmd << " --un #{unal}"
      end
      cmd << " -S /dev/null"

      align = Cmd.new(cmd)
      align.run unal
      hash = {}
      if align.status.success?
        if right
          left[:filtered] = "#{unal}.1"
          right[:filtered] = "#{unal}.2"
          left[:current] = "#{unal}.1"
          right[:current] = "#{unal}.2"
        else
          left[:current] = unal
          left[:filtered] = unal
        end
      else
        msg = "Bowtie2 failed\n#{align.stdout}\n#{align.stderr}"
        raise RuntimeError.new(msg)
      end
    end

    def get_stats
      str=""
      @stats.each do |key, value|
        str << "#{key}\t#{value}\n"
      end
      return str
    end

  end

end
