/* MOLecular DYnamics simulation code, Moldy.
Copyright (C) 1988, 1992, 1993 Keith Refson

This program is free software; you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by the
Free Software Foundation; either version 2, or (at your option) any
later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, 675 Mass Ave, Cambridge, MA 02139, USA.

In other words, you are welcome to use, share and improve this program.
You are forbidden to forbid anyone else to use, share and improve
what you give them.   Help stamp out software-hoarding! */
/*****************************************************************************
 * Parser for transformation expressions a la International Tables (eg       *
 * "x-y,y,3*Z+1/3")							     *
 * BNF form is								     *
 * [trans] :== [const],...,[const]|[sexpr],[sexpr],[sexpr] ;		     *
 * [sexpr] :== [aop][expr]|[expr] ;					     *
 * [expr]  :== [el][aop][expr]|[el] ;					     *
 * [el]    :== [cexpr]*[var]|[cexpr]|[var] ;				     *
 * [aop]   :== +|- ;							     *
 * [cexpr] :== [const]/[const]|[const] ;				     *
 * [var]   :== x|y|z ;							     *
 *  where [const] is built in and detected by lexer.			     *
 *****************************************************************************/

#include <stdio.h>
#include <ctype.h>
#include <stdlib.h>
#include "string.h"
#include <math.h>

#include "sginfo.h"

double 	strtod(const char *, char **);
# define MAXCONST 127
static float constant[MAXCONST];

#define LITERAL 128
#define ISEXP	LITERAL+4
#define IEXPR	LITERAL+7
#define IEL  	LITERAL+10
#define IAOP 	LITERAL+14
#define ICEXP	LITERAL+17
#define IVAR 	LITERAL+20
#define ICONST	LITERAL

#define JTRANS1 1
#define JTRANS2 25
#define JSEXP1	31
#define JSEXP2	34
#define JEXP1	36
#define JEXP2	40
#define JEL1	42
#define JEL2	46
#define JEL3	48
#define JAOP1	50
#define JAOP2	52
#define JCEXP1	54
#define JCEXP2	58
#define JVAR1	60
#define JVAR2	62
#define JVAR3	64

/*
 * Parse tables.  "ps" contains match text and "pstab" master alternatives.
 */

static int ps[] =       {0,ICONST,',',ICONST,',',ICONST,',',ICONST,',',
			ICONST,',',ICONST,',',ICONST,',',ICONST,',',
			ICONST,',',ICONST,',',ICONST,',',ICONST,0,
			ISEXP,',',ISEXP,',',ISEXP,0,
			IAOP,IEXPR,0,IEXPR,0,
			IEL,IAOP,IEXPR,0,IEL,0,
			ICEXP,'*',IVAR,0,ICEXP,0,IVAR,0,
			'+',0,'-',0,
			ICONST,'/',ICONST,0,ICONST,0,
			'x',0,'y',0,'z',0};

int pstab[] =	         {0,JTRANS1,JTRANS2,0,
			 JSEXP1,JSEXP2,0,
			 JEXP1,JEXP2,0,
			 JEL1,JEL2,JEL3,0,
			 JAOP1,JAOP2,0,
			 JCEXP1,JCEXP2,0,
			 JVAR1,JVAR2,JVAR3,0};

#define MATCH   1
#define NOMATCH 0
/*
 * The parser
 */
int parse_trans(char **str, int *pstp, int **parse_rec)
{
   char *strp;
   int  *prec, *pst;

   prec = *parse_rec;
   for(; *pstp; pstp++)			/* Loop over possible alternatives */
   {
      strp = *str;
      *parse_rec=prec;
      *prec = pstp-pstab;
      *(prec+1) = 0;
      for(pst=ps+(*pstp); *pst; pst++)  /* Loop over elements of phrase */
      {
	 if( *pst < LITERAL)
	 {
	    if( *pst != *strp )
	       break;
	    strp++;
	 }
	 else if( *pst == LITERAL )
	 {
	    if( ! (*strp & LITERAL) )
	       break;
	    strp++;
	 }
	 else
	 {
	    (*parse_rec)++;
	    if( ! parse_trans(&strp, pstab+(*pst-LITERAL), parse_rec ))
	       break;
	 }
      }
      if( ! *pst ) 		/* pst==0 if matched whole phrase     */
      {
	 *str = strp;		/* Mark end of text scanned so far    */
	 return MATCH;
      }
   }
   *str = strp;		/* Mark end of text scanned so far    */
   return NOMATCH;
}

/*
 * Lexical analyzer
 */

char *process_line(char *line)
{
   int n = strlen(line);
   char *outline = malloc(n+1);
   char *ip = line, *op = outline;
   int iconst=0;

   if( !outline )
      error("Process_line: Failed to allocate %d bytes of memory\n", n+1);
   
   while( ip < line+n)
   {
      while(*ip == ' ' || *ip == '\t' )
	 ip++;
      
      if(isdigit(*ip) || (*ip == '.' && isdigit(*(ip+1))))
      {
	 constant[iconst++] = strtod(ip, &ip);
	 *op++ = LITERAL;
      }
      else
      {
	 *op = *ip;
	 if( isupper(*op) )
	    *op = tolower(*op);
	 op++; ip++;
      }
   }
   *op = 0;
   return outline;
}

/******************************************************************************
 * make_trans_matrix().  Release parse record and construct transf. matrix.   *
 ******************************************************************************/
int make_trans_matrix(T_RTMx *transformation, int *parse_rec)
{
   int	coeff;
   int		irow = -1, i, j, aflg, iconst = 0;

   for( i=0; i<3*4; i++)
      transformation->a[i] = 0;

   while( *parse_rec )
   {
      switch( pstab[*parse_rec++] )
      {
       case JTRANS1:
	 for(i=0; i<3; i++)
	    for(j=0; j<4; j++)
	       transformation->a[3*i+(j+3)%4] = STBF*constant[iconst++]+0.5;
	 break;
       case JTRANS2:
	 break;
       case JSEXP1:
	 irow++;
	 break;
       case JSEXP2:
	 coeff = STBF;
	 irow++;
	 break;
       case JEXP1:
       case JEXP2:
	 break;
       case JEL1:
       case JEL3:
	 aflg = 0;
	 break;
       case JEL2:
	 aflg = 1;
	 break;
       case JAOP1:
	 coeff = STBF;
	 break;
       case JAOP2:
	 coeff = -STBF;
	 break;
       case JCEXP1:
	 coeff *= constant[iconst++];
	 coeff /= constant[iconst++];
	 if( aflg )
	    transformation->s.T[irow] += coeff;
	 break;
       case JCEXP2:
	 coeff *= (int)(STBF*constant[iconst++]+0.5);
	 coeff /= STBF;
	 if( aflg )
	    transformation->s.T[irow] += coeff;
	 break;
       case JVAR1:
	 transformation->s.R[3*irow+0] += coeff/STBF;
	 break;
       case JVAR2:
	 transformation->s.R[3*irow+1] += coeff/STBF;;
	 break;
       case JVAR3:
	 transformation->s.R[3*irow+2] += coeff/STBF;;
	 break;
       default:
	 error("Unknown entry in parse record, %d\n", pstab[*--parse_rec]);
      }
   }
   return 1;
}

/******************************************************************************
 * transformation_matrix()   Driver function for above			      *
 ******************************************************************************/
int	transformation_matrix(char *buf, T_RTMx *trans_matrix)
{
   int parse_rec[80],*prp;
   char	*tfbuf, *tfp;

   tfbuf = tfp = process_line(buf);
   if( tfbuf == 0 )
      return NOMATCH;
   prp = parse_rec;
   if( parse_trans(&tfp,pstab+1,&prp) && 
      (*tfp == 0 || *tfp == ':' || *tfp == ';') )
   {
      if( ! make_trans_matrix(trans_matrix,parse_rec) )
	 return NOMATCH;
#ifdef DEBUG
      {int i;
      for(i=0; i<3; i++)
	 debug(3,"%8f %8f %8f %8f\n",
		 trans_matrix[i][0], trans_matrix[i][1],
		 trans_matrix[i][2], trans_matrix[i][3]);
      }
#endif
   }
   else
      return NOMATCH;
   (void)free(tfbuf);
   return MATCH;
}

