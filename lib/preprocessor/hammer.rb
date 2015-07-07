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
      dir = "hammer-#{left[:name]}-"
      dir << "#{left[:type]}-" if left[:name]!=left[:type]
      dir << "#{left[:rep]}"
      dir = File.join(@outdir, dir)
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
      @phred = detect_phred(left[:current])
      cmd << "--phred-offset #{@phred} "
      hammer_cmd = Cmd.new(cmd)
      yaml_file = File.join(dir, "corrected", "corrected.yaml")
      unless File.exist?(yaml_file)
        hammer_cmd.run
        if !hammer_cmd.status.success?
          msg = "BayesHammer failed\n#{hammer_cmd.stdout}\n#{hammer_cmd.stderr}"
          raise RuntimeError.new(msg)
        end
      end
      yaml = YAML.load_file(yaml_file)
      cat_cmd = "cat "
      left[:prehammer] = left[:current]
      right[:prehammer] = right[:current]
      left[:current] = yaml[0]["left reads"][0]
      right[:current] = yaml[0]["right reads"][0]

      yaml.each do |item|
        if item.key?("single reads")
          cat_cmd << "#{item["single reads"]} "
        end
      end
      out = "#{dir}/single_reads.fq"
      cat_cmd << " > #{out}"
      unless File.exist?(out)
        cat = Cmd.new(cat_cmd)
        cat.run
      end
      left[:unpaired] = out
      # right[:unpaired] = out
    end

    def run_batch left, right # array of file names
      left_files = []
      right_files = []
      single = []
      left.each do |info|
        left_files << info[:current]
        single << info[:unpaired] if info[:unpaired] and info[:unpaired].length > 0
      end
      right.each do |info|
        right_files << info[:current]
        single << info[:unpaired] if info[:unpaired] and info[:unpaired].length > 0
      end
      gem_dir = Gem.loaded_specs['preprocessor'].full_gem_path
      dir = "hammer"
      dir = File.join(@outdir, dir)
      phred = detect_phred left.first[:current]
      dataset = "[\n  {\n    orientation: \"fr\",\n"
      dataset << "    type: \"paired-end\",\n"
      dataset << "    left reads: [\n"
      dataset << "      #{left_files.join(",\n      ")}\n"
      dataset << "    ],\n    right reads: [\n"
      dataset << "      #{right_files.join(",\n      ")}\n"
      dataset << "    ]\n  }"
      if single and single.length > 0
        dataset << ",\n  {\n"
        dataset << "    type: \"single\",\n"
        dataset << "    single reads: [\n"
        dataset << "      #{single.join(",\n      ")}\n"
        dataset << "    ]\n  }"
      end
      dataset << "\n]\n"
      dataset_file = "#{@outdir}/dataset.yaml"
      File.open(dataset_file, "wb") { |io| io.write dataset }
      puts "wrote dataset file to #{dataset_file}"
      cmd = "#{@spades} "
      cmd << "--dataset #{dataset_file} "
      cmd << "--only-error-correction "
      cmd << "--disable-gzip-output "
      cmd << "-t #{@threads} "
      cmd << "-m #{@memory} "
      cmd << "-o #{dir} "
      cmd << "--phred-offset #{phred} "
      puts cmd
      b = Cmd.new cmd
      yaml_file = "#{dir}/corrected/corrected.yaml"
      b.run yaml_file
      yaml = YAML.load_file(yaml_file)
      left_output = yaml.first["left reads"]
      right_output = yaml.first["right reads"]
      single_output = []
      yaml.each do |hash|
        if hash and hash.key?("single reads")
          hash["single reads"].each do |file|
            single_output << file
          end
        end
      end


      left.zip(left_output) do |a, b|
        a[:prehammer] = a[:current]
        a[:current] = b
        a[:processed][:correction] = "bayeshammer"
      end
      right.zip(right_output) do |a, b|
        a[:prehammer] = a[:current]
        a[:current] = b
        a[:processed][:correction] = "bayeshammer"
      end
      left.zip(single_output).each do |a,b|
        a[:unpaired] = b
      end
    end

    def stats info
      # open :prehammer and :current and compare bases
      # hopefully the reads will be the same length
      # create a list of base position where corrections were made
      @errors = Array.new(100,0)
      @error_qualities = []
      before = File.open(info[:prehammer])
      after = File.open(info[:current])
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
      qe = "#{@outdir}/quality_and_errors.txt"
      if @phred==33
        max = 74
      elsif @phred==64
        max = 105
      end
      min = max - 41
      unless File.exist?(qe)
        File.open(qe, "wb") do |out|
          out.write "base\t#{(min..max).reduce([]) {|list, i| list << i}.join("\t")}\n"
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
