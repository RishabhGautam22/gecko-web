FILENAME=="dico" {code[$1]=$2}
FILENAME=="scheme.old"{
  for(i=1; i<=NF; i++) 
     {
        if (i<NF) 
           { if (substr($i,2,6) in code)
              { printf "%s\t",code[substr($i,2,6)]}
             else
              {printf "%s\t",$i}              
           }
        else
           {printf "%s\n",$i}
            
      }
}
