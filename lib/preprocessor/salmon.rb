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
    end

    def get_stats
      str=""
      return str
    end

  end

end
