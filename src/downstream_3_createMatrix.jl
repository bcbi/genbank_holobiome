#=
This program takes in the names of a phrase and file names mapped host x speices, and output, and creates a matrix of hosts x species associated with phrase
**Double check the column names of CSV match column attributes in code**

=#
using DelimitedFiles, CSV, Suppressor


function hostSpeciesMappings_set(phrase, hostFile, grouping, resolution)

	hostNames=Set(String[])
	speciesNames=Set(String[])
	speciesHost=Dict{String,Array{String}}()
	println("in the hostSpeciesMapping")
	removeSpecies=["Chordata","Arthropoda","Mollusca","fish metagenome","Echinodermata","Cnidaria","Hemichordata","Brachiopoda","activated sludge", "India", "Tea", "unclassified bacterium", "unclassified eukaryote", "unclassified organism", "unclassified circular", "unclassified Antarctic", "unclassified microorganism", "unclassified fungus", "unclassified marine", "unclassified candidate", "unclassified virus", "unclassified phage", "unclassified prokaryote", "unclassified archaeon", "unclassified symbiotic", "pig metagenome", "fish metagenome", "bird metagenome", "gut metagenome", "Northern red-backed", "Duck faeces", "Bovine faeces", "fungal sp.", "Bovine enteric", "termite gut", "mixed culture", "Chicken stool", "Chaerephon bat", "Llama faeces", "lobster gut", "Chamois faeces", "Human feces", "Pacific flying", "Chicken stool-associated", "Platyhelminthes", "unclassified yeast", "Porcine faeces", "Pepper mild", "invertebrate environmental", "Satellite phage", "Deer faeces", "Human gut-associated", "Eidolon bat", "Poales", "Hubei chrysolike", "Rousettus bat", "Kenya bat", "Rhinolophus bat", "eukaryote sp.", "Alphatorquevirus", "Aquabirnavirus", "Annelida", "Nirovirus", "HP38A_virus", "Uroviricota", "Wuhan_large", "Peploviricota"]
	chosenSpecies = ["Crenarchaeota", "Chytridiomycota", "Mucoromycota","Chlorophyta","Streptophyta","Planctomycetes","Elusimicrobia","Candidatus_Saccharibacteria","Tenericutes","Fusobacteria","Cyanobacteria","Dinophyceae","Lentisphaerae","Chlamydiae","Bacteroidetes","Actinobacteria","Acidobacteria","Blastocystis","Verrucomicrobia","Candidatus_Melainabacteria","Chloroflexi","Zoopagomycota","Basidiomycota","Ruminobacillus","Proteobacteria","Synergistetes","Ascomycota","Candidatus_Thermoplasmatota","Aquificae","Firmicutes","Oomycota","Euryarchaeota","Spirochaetes","Fibrobacteres","Deferribacteres", "Microsporidia","Aquabirnavirus", "Uroviricota", "Peploviricota","Artverviricota", "Riboviria_sp.","Anelloviridae","Cossaviricota","Kitrinoviricota","Pisuviricota","Preplasmiviricota","Anellovirus","Bat_fisalivirus","Bat_felisavirus","Cressdnaviricota","Avibirnavirus","Bat_posalivirus","Camula_virus", "Statovirus","Iotatorquevirus","Negarnaviricota","Lenarviricota","Duplornaviricota","Gyrovirus","Thetatorquevirus","Evosea","Heterolobosea","Annelida","Acanthocephala","Phixviricota","Discosea","Fornicata","Nematoda","Cercozoa","Apicomplexa","Tubulinea","Parabasalia","Ciliophora","Euglenozoa","Deinococcus-Thermus","Amoebozoa"]
	open(hostFile) do hostFile
		@suppress begin
			for row in CSV.File(hostFile)
				if grouping=="GROUP"
					groupColumn=row.group
				else
					groupColumn=row.preferred
				end

				if resolution=="GENUS"
					speciesColumn=row.speciesGenus
				elseif resolution=="FAMILY"
					speciesColumn=row.speciesFamily
				elseif resolution=="ORDER"
					speciesColumn=row.speciesOrder
				elseif resolution=="CLASS"
					speciesColumn=row.speciesClass
				elseif resolution=="PHYLUM"
					speciesColumn=row.speciesPhylum
				else
					speciesColumn=row.species
				end
				if (occursin(phrase, groupColumn) & (row.speciesPhylum in chosenSpecies))
					host = cleanHost(row.host_cleaned)
					species = cleanSpeciesName(speciesColumn)
					if !(ismissing(host))
						if (!(species in speciesNames) | !(host in hostNames))
							push!(speciesNames, species)
							push!(hostNames, host)
						end
						if haskey(speciesHost, species) 
							push!(speciesHost[species], host)
						else
							speciesHost[species]=[host]
						end
					end
				end
			end
		end
	end
	return hostNames, speciesNames, speciesHost
end

function cleanSpeciesName(value)
	if occursin(",", value)
		value = replace(value, ","=>"_")
		value = replace(value, ":"=>"_")
		value = replace(value, "["=>"")
		value = replace(value, "]"=>"")
	end
	return value
end

function cleanHost(value)
	if ismissing(value)
		return missing
	end
	if occursin("/",value)
		value = replace(value, "/"=>"_")
	end
	value = replace(value, "."=>"")
	return value
end

function hostSpecieMatrix(hostNames, speciesNames, speciesHost, outputFile, phrase, threshold)
	hostNames=collect(hostNames)
	speciesNames=collect(speciesNames)
	speciesIndices=Dict{String,Array{Int64}}()
	coverageFile = string("hostFiles/", phrase, "_coverage.txt")
	open(coverageFile, "w") do cFile
		for mapping in speciesHost
			specie=String(mapping[1])
			specieIndex=findall(x->x==specie,speciesNames)
			uniqueHosts = unique(mapping[2])
			write(cFile, string(mapping[1], "\n"))
			write(cFile, string(mapping[2], "\n"))
			for host in uniqueHosts
				if haskey(speciesIndices,host)
					speciesIndices[host]=vcat(speciesIndices[host],specieIndex)
				else
					speciesIndices[host]=specieIndex
				end	
			end
		end
	end
	
	updatedHosts=Set(String[])	
	hostFile = string("hostFiles/", phrase, "_hosts.txt")
	open(hostFile, "w") do hFile
		for mapping in speciesIndices
			if length(mapping[2]) >= threshold
				push!(updatedHosts, mapping[1])
				write(hFile,string(mapping[1],"\n"))
			end
		end
	end

	println(phrase," hosts above ",string(threshold*100),"%: ",length(updatedHosts))

	open(outputFile,"w") do outputFile
		write(outputFile,string("#NEXUS\nbegin data;\n"))
		write(outputFile, string("dimensions ntax=", length(updatedHosts), " nchar=", length(speciesNames), ";\n"))
		species=string(speciesNames)
		species=replace(species,", "=>",")
		species=replace(species," "=>"_")
		species=replace(species,","=>" ")
		species=replaceAll(species, ["(", ")","[", "]","{", "}",'"',''', ";","/", ":", "=", "*", "+", "-", "<", ">","@"])
		
		speciesCheck=split(species)
		println(length(speciesCheck))
	
		write(outputFile, string("charlabels ", species, ";\n"))
		write(outputFile, "matrix\n")
		for mapping in speciesIndices
			row = zeros(Int8,length(speciesNames))
			for i in mapping[2]
				row[i]=1
			end
			if mapping[1] in updatedHosts
				row=string(row)
				row =replace(row, "Int8["=>"")
				row =replace(row, ","=>"")
				row =replace(row, "]"=>"")
				host=replace(mapping[1]," "=>"_")
				host=replaceAll(host,["(", ")","[", "]","{", "}",'"',''',",", ";", ":", "=", "*", "+", "-", "<", ">","@"])
				write(outputFile, string(host," ",row,'\n'))
			end
		end
		write(outputFile, ";\n")
		write(outputFile, "end;\n")
	end

end

function replaceAll(str, replacements)
	for r in replacements
		str=replace(str, r=>"")
	end
	return str
end

function main(args)
	if length(args)==6
		hostNames, speciesNames, speciesHost=hostSpeciesMappings_set(uppercase(args[1]),args[2],uppercase(args[4]),uppercase(args[5]))
		hostSpecieMatrix(hostNames, speciesNames, speciesHost, args[3], uppercase(args[1]), parse(Int,args[6]))
		
	else
		println("please enter the names of phrase, and mapped speices to host, and output files as commandline arguments.")
	end

end

main(ARGS)
