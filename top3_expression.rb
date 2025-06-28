# 2024-08-16
# LGT expression file process
if ARGV.size <3
	puts "ruby top3_expression.rb target_gens.txt profile.tsv sample.tsv"
	abort
end

file_target = ARGV.shift
file_profile =ARGV.shift
file_sample =ARGV.shift


# load target gene to count
# target file contains 1 column each line
targets=Hash.new
File.open(file_target) do |tgt|
	tgt.each do |line|
		targets[line.chomp]=1
	end
end

## profile file
##chrom	chromStart	chromEnd	name	score	strand	geneId	geneType	expCount	expScores
#   0     1           2           3       4       5       6        7          8            9
collected=[]
File.open(file_profile) do |prf|
	prf.each do |line|
		cols = line.split(/\t/)
		chr = cols[0]
		start = cols[1]
		stop = cols[2]
		transcript = cols[6]#cols[3]		
		gene = cols[3] #cols[6]
		profile = cols[9]
		strand = cols[5]
		if targets.has_key?(gene)
			collected << [chr,start,stop,strand,transcript,gene,profile]
			# 5-gene
		end
	end
end


# load sample
# tissue file
sample_anno=[]
File.open(file_sample) do |sp|
	sp.each do |line|
		sample_anno<< line.chomp
	end
end

def eval_float(a)
	if a.include?("e")
		t = a.split(/e/)
		return t[0].to_f * (10 ** t[1].to_i)
	end
	return a.to_f
end

def max_3(ee)
	l=ee.chomp.split(/,/)
	max_i1=0
	max_i2=0
	max_i3=0
	l.each_index do |i|
		if eval_float(l[i])>eval_float(l[max_i1])
			max_i1=i
		end
	end

	l.each_index do |i|
		if i==max_i1
			next
		end

		if eval_float(l[i])>eval_float(l[max_i2])
			max_i2=i
		end
	end

	l.each_index do |i|
		if i==max_i1 or i==max_i2
			next
		end

		if eval_float(l[i])>eval_float(l[max_i3])
			max_i3=i
		end
	end
	return [max_i1,max_i2,max_i3]

end
puts "Chromosome\tStart\tStop\tStrand\tID\tGene\tTop1Tissue\tTop2Tissue\tTop3Tissue"
collected.each do |line|
	# [chr,start,stop,strand,transcript,gene,profile]
	chr = line[0]
	start= line[1]
	stop =line [2]
	strand =line[3]
	transcript=line[4]
	gene = line[5]
	profile =line[6]

	l=profile.chomp.split(/,/)
	m1,m2,m3 = max_3(profile)

	o1=sample_anno[m1]
	o2=sample_anno[m2]
	o3=sample_anno[m3]

	if l[m2]=="0"
		o2 = " "
	end
	if l[m3]=="0"
		o3 = " "
	end
	if targets.has_key?(gene)
		puts "#{chr}\t#{start}\t#{stop}\t#{strand}\t#{transcript}\t#{gene}\t#{o1}\t#{o2}\t#{o3}"
	end
end

