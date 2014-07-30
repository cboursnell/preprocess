require 'bindeps'

module Preprocessor

  class Skewer

    def initialize(outdir, threads, end_quality, mean_quality, min_length)
      @outdir = outdir
      @end_quality = end_quality
      @mean_quality = mean_quality
      @min_length = min_length
      @filter = true
      gem_dir = Gem.loaded_specs['preprocessor'].full_gem_path
      gem_deps = File.join(gem_dir, 'deps', 'skewer.yaml')
      Bindeps.require gem_deps
    end

    def run left, right=nil
      if right
        trim_cmd = "skewer-0.1.117-linux-x86_64"
        trim_cmd << " -m pe"
        trim_cmd << " -q #{@end_quality}"
        trim_cmd << " -Q #{@mean_quality}"
        trim_cmd << " -l #{@min_length}"
        trim_cmd << " -n " if @filter
        trim_cmd << " -o #{@outdir}/#{left[:name]}-#{left[:type]}-#{left[:rep]}"
        trim_cmd << " #{left[:current]}"
        trim_cmd << " #{right[:current]}"
        cmd = Cmd.new(trim_cmd)
        cmd.run
        out = cmd.stdout
        if out =~ /log has been saved to \"(.*)\"/
          log = $1
          File.open("#{log}").each_line do |line|
            if line =~ /^trimmed:\s*(.+)$/
              filenames = $1
              files = filenames.split(/,\s*/)
              left[:current] = files[0]
              right[:current] = files[1]
            end
          end
        end
      else
      end
    end

  end

end
