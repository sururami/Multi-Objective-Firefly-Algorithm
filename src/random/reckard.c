/* 
   ForwardRandomUniform()
   Skip ahead in the random sequence.   
   by Stan Reckard  
*/
void ForwardRandomUniform(long forward)
{
   double uni;

   /* Make sure the initialisation routine has been called */
   if (!test)
      RandomInitialise(1802,9373);

   while (forward--) {
       uni = u[i97-1] - u[j97-1];
       if (uni <= 0.0)
          uni++;
       u[i97-1] = uni;
       i97--;
       if (i97 == 0)
          i97 = 97;
       j97--;
       if (j97 == 0)
          j97 = 97;
       c -= cd;
       if (c < 0.0)
          c += cm;
   }
}

/* 
   BackupRandomUniform()
   Backup in the random sequence 'backup' times.   
   by Stan Reckard  
*/
void BackupRandomUniform(long backup) 
{
   double uni, uniAlt, prev;
   double cTmp;

   while (backup--) {
      if (c >= cm)
          c -= cm;

      c += cd;

      if (j97 == 97)
         j97 = 0;
      j97++;

      if (i97 == 97)
         i97 = 0;
      i97++;

      uni = u[i97-1];
      uniAlt = uni - 1;

      prev = uni + u[j97-1];
      if ((prev > 0.0F) && (prev < 1.0F))
         u[i97-1] = prev;
      else
         u[i97-1] = uniAlt + u[j97-1];;
   }
}

/* 
   UnRandomUniform()
   Backup in the random sequence.   
   by Stan Reckard 
*/
double UnRandomUniform(void) 
{
   double uni, uniAlt, prev;
   double cTmp;

   if (c >= cm)
      c -= cm;

   c += cd;

   if (j97 == 97)
      j97 = 0;
   j97++;

   if (i97 == 97)
      i97 = 0;
   i97++;

   uni = u[i97-1];
   uniAlt = uni - 1;

   prev = uni + u[j97-1];
   if ((prev > 0.0F) && (prev < 1.0F))
      u[i97-1] = prev;
   else {
      u[i97-1] = uniAlt + u[j97-1];
   }
   /* RandomUniform() has been completely undone at this point.  */

   /*
      Now get the random# that was last retrieved to 
		return without altering the random sequence.
      uni holds old value of u[i97-1]
   */
   cTmp = c;
   cTmp -= cd;
   if (cTmp < 0.0F)
      cTmp += cm;
   uni -= cTmp;
   if (uni < 0.0F)  uni++;

   return uni;  /* prev random# returned */
}


/* rand24()    24-bit precision   */
unsigned int rand24(void) 
{ 
   return (unsigned int)(RandomUniform() * 4096 * 4096); 
}

/* unRand24()    24-bit precision   */
unsigned int unRand24(void) 
{ 
   return (unsigned int)(UnRandomUniform() * 4096 * 4096); 
}

