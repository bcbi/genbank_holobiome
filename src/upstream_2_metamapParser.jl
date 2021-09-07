#=
This program takes in metamap output and creates a two CSV formated in the following:
CSV 1: Original String, MetaMap's Output, Semantic Type
CSV 2: Original string, CUI, MetaMap Value, MetaMap's assigned Semantic Type
=#

using DelimitedFiles


function writeToCSV(parserOutputFile, cuiOutputFile, phrase, mappings, semanticTypes)
	open(parserOutputFile, "a") do parserOutputFile
		open(cuiOutputFile, "a") do cuiOutputFile
			for m in mappings
				mappingSemanticTypes = split(m[2],'|')
				for s in mappingSemanticTypes
					if s in semanticTypes

						preferred=split(m[1],'|')
						fullPreferred=""
						for i=1:length(preferred)
							concept=split(preferred[i], ':')
							fullPreferred=string(fullPreferred," ", concept[2])
							writedlm(cuiOutputFile, [phrase concept[1] concept[2] mappingSemanticTypes[i]], ',')
						end

						fullPreferred=fullPreferred[2:end]
						writedlm(parserOutputFile, [phrase fullPreferred m[2]], ',')

						break
					end
				end 
			end
		end
	end

end 

function parseMetamapFile(MetaMapOutputFile, semanticTypeFile, parserOutputFile, cuiOutputFile)
	if !(isfile(parserOutputFile))
		open(parserOutputFile, "a") do parserOutputFile
			writedlm(parserOutputFile, ["original" "preferred" "semanticType"], ',')
		end
	end

	if !(isfile(cuiOutputFile))
		open(cuiOutputFile, "a") do cuiOutputFile
			writedlm(cuiOutputFile, ["original" "CUI" "concept" "semanticType"], ',')
		end
	end
	semanticTypes = readdlm(semanticTypeFile)  
	open(MetaMapOutputFile) do MetaMapOutputFile
		phrase=""
		prevMapping=""
		mappings=Dict{String,String}()
		prevLnMetaMap=false
		sameMapping=false
		lines = readlines(MetaMapOutputFile)
		numLines=length(lines)
		firstMapping=false
			
		for i in 1:numLines
			ln=lines[i]
			if occursin("Phrase", ln)	
				phraseArray=split(ln,' ')
				phraseArray=phraseArray[2:end]
				phrase=join(phraseArray,' ')
				if occursin(',',ln)
					phrase=phrase[1:end-1]
				end
				firstMapping=true
			elseif occursin("Meta Mapping", ln) 
				prevLnMetaMap=true
				sameMapping=false
			elseif prevLnMetaMap 
				prevMapping=match(r"C.*[\[]",ln).match[1:end-2]
				semanticType=match(r"(?<=\[).+?(?=\])",ln).match
				mappings[prevMapping]=semanticType
				prevLnMetaMap=false
				sameMapping=true
				if i==numLines
					writeToCSV(parserOutputFile, cuiOutputFile,  phrase, mappings, semanticTypes)
				end
			elseif sameMapping & (ln!="") & !(occursin("Processing", ln))
				mapping=match(r"C.*[\[]",ln).match[1:end-2]
				semanticType=match(r"(?<=\[).+?(?=\])",ln).match
				prevSemanticType=mappings[prevMapping]			
				delete!(mappings,prevMapping)
				newMatch=string(prevMapping,"|",mapping)
				newSemanticType=string(prevSemanticType,"|", semanticType)
				mappings[newMatch]=newSemanticType
				prevMapping=newMatch
				if i==numLines
					writeToCSV(parserOutputFile, cuiOutputFile, phrase, mappings, semanticTypes)
				end
			elseif (ln=="") & (phrase!="")
				if mappings!=Dict{String,String}()
					writeToCSV(parserOutputFile, cuiOutputFile, phrase, mappings, semanticTypes)	
				end
				phrase=""
				mappings=Dict{String,String}()
				prevLnMetaMap=false
				sameMapping=false			
			end
		end
	end

end

function main(args)
   
	if length(args)==4
	    parseMetamapFile(args[1], args[2], args[3], args[4])     
	else
		println("please enter the names of metamap output, semantic type files, parsed metmap mapping output, cui mappings as commandline arguments.")
	end
	
	
end

main(ARGS)