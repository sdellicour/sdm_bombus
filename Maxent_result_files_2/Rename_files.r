files1 = list.files()
files2 = list.files()

files1 = files1[which(!grepl("_2.",files1))]
files2 = files2[which(!grepl("_2.",files2))]

for (i in 1:length(files1))
	{
		file1 = files1[i]
		# file2 = gsub("\\.RData","\\_3.RData",file1)
		file2 = gsub("\\.csv","\\_3.csv",file1)
		file2 = gsub("\\.pdf","\\_3.pdf",file2)
		file2 = gsub("\\.asc","\\_3.asc",file2)
		files2[i] = file2
		file.rename(file1,file2)
	}
