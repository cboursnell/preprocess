require 'bindeps'
require 'csv'
require 'fixwhich'

module Preprocessor

  class Salmon

    def initialize(outdir, threads, reference)
      @outdir = outdir
      @threads = threads
      @reference = File.expand_path(reference)

      gem_dir = Gem.loaded_specs['preprocessor'].full_gem_path
      gem_deps = File.join(gem_dir, 'deps', 'salmon.yaml')
      Bindeps.require gem_deps
      @salmon = Which.which('salmon').first
    end

    def run left, right=nil
      if left[:sam]
        sam = left[:sam]
        salmon_out = File.join(@outdir, "salmon_#{File.basename(left[:file])}")

        cmd = "#{@salmon} quant"
        cmd << " --libType IU"
        cmd << " --alignments #{sam}"
        cmd << " --targets #{@reference}"
        cmd << " --threads #{@threads}"
        cmd << " --useErrorModel"
        cmd << " --output #{salmon_out}"

        salmon = Cmd.new(cmd)
        salmon.run
        if !salmon.status.success?
          msg = "Something went wrong with salmon\n#{salmon.stdout}\n"
          msg << "#{salmon.stderr}"
          raise RuntimeError.new(msg)
        end
        left[:salmon] = File.join(salmon_out, "quant.sf")
      elsif left[:current] =~ /fastq/ or left[:current] =~ /fq/
        salmon_out = File.join(@outdir, "salmon_#{File.basename(left[:file])}")
        index = "#{File.basename(@reference, File.extname(@reference))}_index"
        index = File.join(@outdir, index)
        index_cmd = "#{@salmon} index"
        index_cmd << " --index #{index}" # create index
        index_cmd << " --transcripts #{@reference}"
        index_cmd << " --threads #{@threads}"
        make_index = Cmd.new(index_cmd)
        if !File.exist?(File.join(index, "bwaidx.bwt"))
          make_index.run
        end

        cmd = "#{@salmon} quant"
        if left and right
          cmd << " --libType IU"
          cmd << " -1 #{left[:current]}"
          cmd << " -2 #{right[:current]}"
        else
          cmd << " --libType U"
          cmd << " -r #{left[:current]}"
        end
        cmd << " --index #{index}"
        cmd << " --output #{salmon_out}"
        salmon = Cmd.new(cmd)
        salmon.run
        if !salmon.status.success?
          msg = "Something went wrong with salmon\n#{salmon.stdout}\n"
          msg << "#{salmon.stderr}"
          raise RuntimeError.new(msg)
        end
        left[:salmon] = File.join(salmon_out, "quant.sf")

      end

    end

    def get_stats
      str=""
      return str
    end

  end

end
