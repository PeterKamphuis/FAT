;+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
;+
;
;*NAME:	SUM
;
;*CLASS:
;
;*CATEGORY:
;
;*PURPOSE:	Total up an array over one of its dimensions.
;
;*CALLING SEQUENCE:
;	RESULT = SUM(ARRAY,DIMENSION)
;
;*PARAMETERS:
; INPUTS:
;	ARRAY = Input array.  May be any type except string.
;	DIMENSION = Dimension to do total over.
; OUTPUTS:
;	The result is an array with all the dimensions of the input array
;	except for the dimension specified, each element of which is the total
;	of the corresponding vector in the input array. 
;
;*SYSTEM VARIABLES USED:
;
;*COMMON BLOCKS:
;
;*INTERACTIVE INPUT:
;
;*SUBROUTINES CALLED:
;
;*FILES USED:
;
;*SIDE EFFECTS:
;
;*EXAMPLES:
;	For example, if A is an array with dimensions of (3,4,5), then the
;	command B = SUM(A,1) is equivalent to
;
;			B = FLTARR(3,5)
;			FOR J = 0,4 DO BEGIN
;				FOR I = 0,2 DO BEGIN
;					B(I,J) = TOTAL( A(I,*,J) )
;				ENDFOR
;			ENDFOR
;
;*RESTRICTIONS:
;	Dimension specified must be valid for the array passed; otherwise the
;	input array is returned as the output array.
;
;*NOTES:
;
;*PROCEDURE:
;
;*MODIFICATION HISTORY:
;	William Thompson	Applied Research Corporation
;	July, 1986		8201 Corporate Drive
;				Landover, MD  20785
;	JKF/ACC 12/31/90  - copied for use at the GHRS DAF
;-
;-------------------------------------------------------------------------------
FUNCTION SUM,ARRAY,DIMENSION
;
IF (N_PARAMS(0) EQ 0) THEN BEGIN
	PRINT,' CALLING SEQUENCE: RESULT= SUM( ARRAY, DIMENSION)'
        retall
end

IF (N_PARAMS(0) LT 2) THEN BEGIN
	PRINT,'*** Function SUM must be called with two parameters:'
	PRINT,'                   ARRAY , DIMENSION'
	RETURN,ARRAY
ENDIF
;
S = SIZE(ARRAY)
N_DIM = S(0)
IF N_DIM EQ 0 THEN $
	MESSAGE,'*** Variable must be an array, name= ARRAY' $
ELSE IF (DIMENSION GE N_DIM) OR (DIMENSION LT 0) THEN   $
	MESSAGE,'*** Dimension out of range, name= ARRAY'   $
ELSE IF N_DIM EQ 1 THEN BEGIN	;Trivial case, equivalent to TOTAL.
	F = TOTAL(ARRAY)
	RETURN,F
ENDIF
;
NEL=N_ELEMENTS(S)
S2 = S( WHERE( INDGEN(NEL) NE FIX(DIMENSION+1) ))		;Set up defining array
S2 = S2(1:N_DIM-1)					;for output variable.
F = FLTARR(S2(0))
FOR IR=1,N_ELEMENTS(S2)-1 DO F = F#FLTARR(S2(IR))
;
TYPE = S(S(0) + 1)					;Chose data type.
CASE TYPE OF						;(At least REAL*4)
	5: F = DOUBLE(F)
	6: F = COMPLEX(F)
	ELSE:
ENDCASE
;
;  Calculate product of sizes of dimensions lower than, equal to, and higher
;  than DIMENSION (NI,NJ,NK respectively).
;
NI = 1
IF DIMENSION GT 0 THEN FOR M = 1,DIMENSION DO NI = NI * S(M)
NJ = S(DIMENSION+1)
NK = 1
IF DIMENSION LT N_DIM-1 THEN FOR M = DIMENSION+2,N_DIM DO NK = NK * S(M)
;
;  Set up index arrays.
;
XIK = LINDGEN( ( NI * NK ) )
XJ = LINDGEN( ( NJ ) )
NIJ = NI*NJ
;
;  Choose whether it is more efficient to loop over NI and NK ...
;
IF NI*NK LT NJ THEN BEGIN
	FOR I = 0,NI-1 DO FOR K = 0,NK-1 DO $
		F(I+NI*K) = TOTAL( ARRAY(I + NI*XJ + NIJ*K) )
;
;  ... or over NJ.
;
END ELSE BEGIN
	XI = XIK MOD NI
	XK = LONG( XIK / NI )
	FOR J = 0,NJ-1 DO F = F + ARRAY(XI + NI*J + NIJ*XK)
ENDELSE
RETURN,F
END
