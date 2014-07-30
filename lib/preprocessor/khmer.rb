require 'bindeps'
require 'which'
include Which

module Preprocessor

  class Khmer

    def initialize(outdir, threads, memory, kmer, cutoff, tables)
      @outdir = outdir
      @threads = threads
      @memory = memory.to_f
      @kmer = kmer
      @cutoff = cutoff
      @tables = tables
      puts @memory
      puts @tables
      @min_table_size = ((@memory/@tables)*1e9).to_i
      puts @min_table_size
      @khmer = which("normalize-by-median.py").first
      if @khmer.empty?
        raise RuntimeError.new("khmer not installed")
      end
      @hashes = []
    end

    def normalize
      paired_file_list = ""
      single_file_list = ""
      paired_name = "paired"
      single_name = "single"
      @hashes.each do |pair|
        hash = pair[0]
        paired_file_list << "#{hash[:current]} "
        single_file_list << "#{pair[0][:unpaired]} " if pair[0][:unpaired]
        single_file_list << "#{pair[1][:unpaired]} " if pair[1][:unpaired]
        paired_name = "#{hash[:name]}.khmered.fq"
        single_name = "#{hash[:name]}.unpaired.fq"
      end
      @paired_outfile = File.join(@outdir, paired_name)
      if single_file_list.length > 1
        @single_outfile = File.join(@outdir, single_name)
      else
        @single_outfile = nil
      end
      # paired
      pair = " -p " if @paired
      cmd = "#{@khmer} #{pair} --quiet "
      cmd << " --ksize #{@kmer} --cutoff #{@cutoff} "
      cmd << " --n_tables #{@tables}  --min-tablesize #{@min_table_size} "
      cmd << " --fault-tolerant --out #{@paired_outfile} "
      cmd << "#{paired_file_list}"
      puts cmd
      khmer_cmd = Cmd.new(cmd)
      khmer_cmd.run
      # unpaired
      if single_file_list.length > 1
        cmd = "#{@khmer} --quiet "
        cmd << " --ksize #{@kmer} --cutoff #{@cutoff} "
        cmd << " --n_tables #{@tables}  --min-tablesize #{@min_table_size} "
        cmd << " --fault-tolerant --out #{@single_outfile} "
        cmd << "#{single_file_list}"
        puts cmd
        khmer_cmd = Cmd.new(cmd)
        khmer_cmd.run
      end
      left_outfile = File.join(@outdir, "#{paired_name}-left.fq")
      right_outfile = File.join(@outdir, "#{paired_name}-right.fq")
      self.deinterleave(@paired_outfile, left_outfile, right_outfile)
      new_data = []
      new_data << {:current => left_outfile, :unpaired => @single_outfile,
                   :name => @hashes[0][0][:name],
                   :type => @hashes[0][0][:type], :rep => 1, :pair => 1}
      new_data << {:current => right_outfile, :unpaired => nil,
                   :name => @hashes[0][0][:name],
                   :type => @hashes[0][0][:type], :rep => 1, :pair => 2}
      return new_data
    end

    def interleave left, right=nil
      if right
        outfile = File.join(@outdir,
                          "#{left[:name]}-#{left[:type]}-#{left[:rep]}.in.fq")
        cmd = "paste #{left[:current]} #{right[:current]} | paste - - - - | "
        cmd << "awk -v FS=\"\\t\" -v OFS=\"\\n\" \'{print(\"@read\"NR\":1\","
        cmd << "$3,$5,$7,\"@read\"NR\":2\",$4,$6,$8)}\' > "
        cmd << " #{outfile}"
        interleaver = Cmd.new(cmd)
        interleaver.run
        left[:current] = outfile
        right[:current] = outfile
        @hashes << [left, right]
        @paired = true
      else
        # don't need to interleave if not paired reads
        @hashes << [left, nil]
        @paired = false
      end
    end

    def deinterleave(file, output_left, output_right)
      #puts "deinterleaving #{file} to make #{output_left} and #{output_right}"
      raise RuntimeError.new("Can't find #{file}") if !File.exist?(file)
      fastq = File.open(file)
      left = File.open("#{output_left}", "w")
      right= File.open("#{output_right}", "w")

      name1 = fastq.readline rescue nil
      seq1 = fastq.readline rescue nil
      plus1 = fastq.readline rescue nil
      qual1 = fastq.readline rescue nil
      name2 = fastq.readline rescue nil
      seq2 = fastq.readline rescue nil
      plus2 = fastq.readline rescue nil
      qual2 = fastq.readline rescue nil

      while name1 != nil and name2 != nil
        left.write(name1)
        left.write(seq1)
        left.write(plus1)
        left.write(qual1)
        right.write(name2)
        right.write(seq2)
        right.write(plus2)
        right.write(qual2)

        name1 = fastq.readline rescue nil
        seq1 = fastq.readline rescue nil
        plus1 = fastq.readline rescue nil
        qual1 = fastq.readline rescue nil
        name2 = fastq.readline rescue nil
        seq2 = fastq.readline rescue nil
        plus2 = fastq.readline rescue nil
        qual2 = fastq.readline rescue nil
      end
      fastq.close
      left.close
      right.close
    end


  end

end
