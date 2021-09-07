#= 
This is an program that can be used to create a text document of input words
for MetaMap from a csv that contains other data. 
**Double check the column names of CSV match column attributes in code**
=#

using DelimitedFiles, CSV, Suppressor

function parseCSV(inputFile, outputFile)
	
	inputValues=String[]

	open(inputFile) do inputFile
		open(outputFile,"a") do outputFile
			@suppress begin
				for row in CSV.File(inputFile)
					if !ismissing(row.value)
						if !(row.value in inputValues)
							push!(inputValues,row.value)
							write(outputFile, string(row.value,",\n"))
						end
					end
				end
			end
		end
	end
	

end


function main(args)
   
	if length(args)==2
	    parseCSV(args[1], args[2])
	else
		println("please enter the names of input and output files as commandline arguments.")
	end
	
	
end

main(ARGS)
