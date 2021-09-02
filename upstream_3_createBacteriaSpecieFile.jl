#=
This file takes in Genbanks csv of of certain source site with locus and csv of species name with locus and creates csv of locus, soure name, and name of species name
(csv without preferred name ie. speciesTissue.csv)
=#


using DelimitedFiles, CSV, Suppressor

function parser(siteFile, speciesFile, outputFile)
	locusSiteDict=Dict{String,String}()
	open(siteFile) do siteFile
		@suppress begin
			for row in CSV.File(siteFile)
				if row.locus!="locus"
					locusSiteDict[row.locus]=row.value
				end
			end
		end
	end

	open(speciesFile) do speciesFile
		open(outputFile, "a") do outputFile
			@suppress begin
				for row in CSV.File(speciesFile)
					if row.locus!="locus"
						if haskey(locusSiteDict, row.locus)
							writedlm(outputFile, [row.locus locusSiteDict[row.locus] row.value], ',')
						end
					end
				end
			end
		end
	end
end


function main(args)

	if length(args)==3
		if !(isfile(args[3]))
			open(args[3], "a") do outputFile
				writedlm(outputFile, ["locus" "original" "species"], ',')
			end
		end

		parser(args[1], args[2], args[3])
	else
		println("please enter the names of site name, species file, and output")
	end

end

main(ARGS)
