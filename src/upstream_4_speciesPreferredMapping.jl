#=
This file creates a csv of metamaps's preffered names of sources
mapped to associated bacteria species. 
(csv with preferred names of sources ie speicesMetaMapTissue.csv)
**Double check the column names of CSV match column attributes in code**
=#

using DelimitedFiles, CSV, Suppressor


function mappingPreferred_to_Species(preferredFile)
	originalToPreferred=Dict{String,Array{String}}()
	open(preferredFile) do preferredFile
		original=""
		preferred=""
		@suppress begin
			for row in CSV.File(preferredFile)
				original=row.original
				preferred=row.preferred
				
				if !haskey(originalToPreferred, original)
					originalToPreferred[original]=[preferred]
				else
					if !(preferred in originalToPreferred[original])
						push!(originalToPreferred[original],preferred)
					end
				end
				
			end
		end
	end
	return originalToPreferred
end


function speciesMapping(originalToPreferred,speciesFile, outputFile)
	open(speciesFile) do speciesFile
		open(outputFile, "a") do outputFile
			@suppress begin
				for row in CSV.File(speciesFile)
					if haskey(originalToPreferred, row.original)
						preferredMappings=originalToPreferred[row.original]
						for mapping in preferredMappings
							writedlm(outputFile, [row.locus row.original mapping row.species], ',')
						end
					end
				end
			end
		end

	end

end

function main(args)
   
	if length(args)==3
		open(args[3], "a") do outputFile
			writedlm(outputFile, ["locus" "original" "preferred" "species"], ',')
		end
	    originalToPreferred=mappingPreferred_to_Species(args[1]) 
	    speciesMapping(originalToPreferred,args[2], args[3])
	    
	else
		println("please enter the names of preferred source names, species, and output files as commandline arguments.")
	end
	
	
end

main(ARGS)
