require 'yaml'

module Preprocessor

  class Hammer

    def initialize(outdir, threads, memory)
      @outdir = outdir
      @threads = threads
      @memory = memory

      gem_dir = Gem.loaded_specs['preprocessor'].full_gem_path
      @spades = File.join(gem_dir, 'bin', 'SPAdes-3.5.0-Linux', 'bin', 'spades.py')
      if !File.exist?(@spades)
        url = "http://spades.bioinf.spbau.ru/release3.5.0/SPAdes-3.5.0-Linux.tar.gz"
        wget_cmd = "wget #{url} -O #{File.join(gem_dir,"bin","spades.tar.gz")}"

        archive = File.join(gem_dir,"bin","spades.tar.gz")
        tar_cmd = "tar xzf #{archive}"
        tar_cmd << " -C #{File.join(gem_dir,"bin")}"
        download = Cmd.new(wget_cmd)
        download.run
        if download.status.success?
          unpack = Cmd.new(tar_cmd)
          unpack.run
          if !unpack.status.success?
            msg = "unpacking spades.tar.gz failed\n#{unpack.stderr}"
            raise RuntimeError.new msg
          end
          File.delete(archive)
        else
          msg = "download of spades.tar.gz failed\n"
          msg << download.stdout
          msg << download.stderr
          raise RuntimeError.new msg
        end
      end

    end

    def run left, right=nil
      gem_dir = Gem.loaded_specs['preprocessor'].full_gem_path
      dir = File.join(@outdir,
                      "hammer-#{left[:name]}-#{left[:type]}-#{left[:rep]}")
      cmd = "#{@spades} "
      cmd << "-1 #{left[:current]} "
      cmd << "-2 #{right[:current]} "
      cmd << "--pe1-s #{left[:unpaired]} " if left[:unpaired]
      cmd << "--pe2-s #{right[:unpaired]} " if right[:unpaired]
      cmd << "--only-error-correction "
      cmd << "--disable-gzip-output "
      cmd << "-t #{@threads} "
      cmd << "-m #{@memory} "
      cmd << "-o #{dir} "
      cmd << "--phred-offset #{detect_phred(left[:current])} "
      puts cmd
      hammer_cmd = Cmd.new(cmd)
      unless Dir.exist? dir
        hammer_cmd.run
        if !hammer_cmd.status.success?
          msg = "BayesHammer failed\n#{hammer_cmd.stdout}\n#{hammer_cmd.stderr}"
          raise RuntimeError.new(msg)
        end
      end
      yaml = YAML.load_file("#{dir}/corrected/corrected.yaml")
      cat_cmd = "cat "
      left[:prehammer] = left[:current]
      right[:prehammer] = right[:current]
      left[:current] = yaml[0]["left reads"][0]
      right[:current] = yaml[0]["right reads"][0]

      if yaml[0] and yaml[0]["single reads"]
        yaml[0]["single reads"].each do |single|
          cat_cmd << " #{single} "
        end
        out = "#{dir}/left_unpaired.fq"
        cat_cmd << " >  #{out}"
        cat = Cmd.new(cat_cmd)
        cat.run
        left[:unpaired] = "#{out}"
      end
      if yaml[1] and yaml[1]["single reads"]
        yaml[1]["single reads"].each do |single|
          cat_cmd << " #{single} "
        end
        out = "#{dir}/right_unpaired.fq"
        cat_cmd << " >  #{out}"
        cat = Cmd.new(cat_cmd)
        cat.run
        right[:unpaired] = "#{out}"
      end

    end

    def stats file
      # open :prehammer and :current and compare bases
      # hopefully the reads will be the same length
      # create a list of base position where corrections were made
      @errors = Array.new(100,0)
      @error_qualities = []
      before = File.open(file[:prehammer])
      after = File.open(file[:current])
      name1 = before.readline
      name2 = after.readline
      while name1 and name2
        seq1 = before.readline
        seq2 = after.readline
        plus1 = before.readline
        plus2 = after.readline
        qual1 = before.readline
        qual2 = after.readline

        while name1 and name1 != name2
          # :current is always going to be a subset of :prehammer
          # so just scan forwards in :prehammer until we find
          # a read with a matching name
          name1 = before.readline rescue nil
          seq1 = before.readline rescue nil
          plus1 = before.readline rescue nil
          qual1 = before.readline rescue nil
        end
        if seq1 and seq1.length == seq2.length
          (0..seq1.length-1).each do |i|
            if seq1[i]!=seq2[i]
              @errors[i] ||= 0
              @errors[i] += 1
              @error_qualities[i] ||= []
              @error_qualities[i][qual1[i].ord] ||= 0
              @error_qualities[i][qual1[i].ord] += 1
            end
          end
        end

        name1 = before.readline rescue nil
        name2 = after.readline rescue nil
      end
    end

    def get_stats
      File.open("#{@outdir}/quality_and_errors.txt", "wb") do |out|
        out.write "base\t#{(32..105).reduce([]) {|list, i| list << i}.join("\t")}\n"
        @error_qualities.each_with_index do |list, index|
          unless list.nil?
            row = []
            row << index
            (32..105).each do |i|
              list[i].nil? ? row << 0 : row << list[i]
            end
            out.write "#{row.join("\t")}\n"
          end
        end
      end
      tot = @errors.reduce(0) { |sum, c| sum += c }
      str = ""
      @errors.each_with_index do |count, length|
        str << "#{length}\t#{count}\t#{(100 * count / tot.to_f).round(2)}%\n"
      end
      str
    end

    def detect_phred filename
      file = File.open(filename)
      c = 0
      scores={}
      while c < 4000
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
