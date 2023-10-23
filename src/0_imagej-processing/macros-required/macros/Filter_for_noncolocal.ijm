//This is to filter for non-colocalised molecules
//Table must be called "Results"
for (i = 0; i < nResults; i++) 
   {
   	d = getResult("distance", i);
if (d>0) {
   	IJ.deleteRows(i,i);
   	i = 0;
   	}
   }

