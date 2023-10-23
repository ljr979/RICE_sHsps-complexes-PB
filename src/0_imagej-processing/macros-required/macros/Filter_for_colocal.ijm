//This is to filter for colocalised molecules
//Table must be called "Results"
for (i = 0; i < nResults; i++) 
   {
   	d = getResult("distance", i);
if (d==-1) {
   	IJ.deleteRows(i,i);
   	i = 0;
   	}
   }

