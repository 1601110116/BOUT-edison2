FUNCTION array_to_gkvsd, array, _Extra=Extra



    ndims = (SIZE(array))[0]
            
    
CASE ndims OF 

  1  :  objStructure = {GKVs1D}
  2  :  objStructure = {GKVs2D}
  3  :  objStructure = {GKVs3D}
  4  :  objStructure = {GKVs4D}
ENDCASE

    objStructure.FileID = ''
    objStructure.RunID = ''
    objStructure.CodeName = ''   
    objStructure.CodePI = ''
    
    objStructure.mnemonic = 'set mnemonic'
    objStructure.Title = 'set title'
    objStructure.Indices = PTR_NEW(REPLICATE('*', ndims))   
    
; 
; make grids
;    




  gridStructure = {Grid}                
  gridStructure.Mnemonic  = 'Grid' 
  gridStructure.Title = 'Grid'
   
  
FOR i = 0,ndims-1 DO BEGIN
  gridStructure = {Grid}                
  gridStructure.Mnemonic  = 'Grid' 
  gridStructure.Title = 'Grid'
  gridStructure.values = PTR_NEW(INDGEN((SIZE(array))[i+1])) 
  gridStructure.Range = [0., (SIZE(array))[i+1]-1]
  gridStructure.irange  = [0, (SIZE(array))[i+1]-1]  
  i2 = i + 1
  command = 'objStructure.grid' + STRTRIM(i2,1) + ' = GKVsd_GridCopy(gridStructure)'  
  OK = EXECUTE(command)  
  PTR_FREE, gridStructure.values
ENDFOR

 

    objStructure.values = PTR_NEW(array)
CASE ndims OF 
  1 : outputStructure = OBJ_NEW('GKVs1D', objStructure)
  2 : outputStructure = OBJ_NEW('GKVs2D', objStructure)
  3 : outputStructure = OBJ_NEW('GKVs3D', objStructure)
  4 : outputStructure = OBJ_NEW('GKVs4D', objStructure)
ENDCASE

RETURN, outputstructure
END