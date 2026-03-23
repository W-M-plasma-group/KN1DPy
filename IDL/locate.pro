;+	LOCATE.PRO
;NAME:		LOCATE
;PURPOSE:	Do a binary search for the index(es) in a sorted table.
;	Finds Index such that Value is between Table(Index) and Table(Index+1).
;	If Table ascends then:
;		Table(Index) LE Value LT Table(Index+1)
;		Index is -1 only for Value LT Table(0) -- too low
;		Index is n-1 only for Value GE Table(n-1) -- high or too high
;	For a descending series, result is (n-1) - the_ascending_series.
;CATEGORY:	Searching and sorting.
;CALL:		Index = LOCATE(Table, Value)
;INPUTS:	Table		sorted list of number/text values (a vector)
;		Value		scalar/array of number/text to look up
;OPTIONAL INPUTS:
;KEYWORDS:
;OUTPUTS:	Index		the index(es) of Value in Table, from 0 to n.
;				Long integer of same shape as Value.
;OPTIONAL OUTPUTS:
;RESTRICTIONS:	Must compare text with text, numbers with numbers
;PROCEDURE:	Adapted and vectorized from Numerical Recipes.
;		Approx. worst timing: (Ceiling(Log2(N_ELEMENTS*(Table))+1)
;		* vector(divide, 2*add, 4*compare, 2*WHERE) on Value
;EXAMPLES:	j = LOCATE([1,2,3,4], 3.5) returns 2L.
;		j = LOCATE([4,3,2,1], [.5,1,4,4.5]) is [3L,3L,0L,-1]
;NOTE:		Should be useful in constructing reverse lookup
;		like arbitrary axis subscripting.
;HISTORY:	09-Apr-1993	Ken Klare, LANL P-4 (c)1993
;		5-Apr-1993	KAK	use MAKE_ARRAY
;	7-Sep-1993	KAK logic for descending. (Labombard@edge1.pfc.mit)
;	7-Sep-1993	KAK allow scalar lookup
;-
FUNCTION LOCATE, Table, Value

 sz = SIZE(Value)					;output shape
 n = N_ELEMENTS(Table)					;size of the table
 asc = Table(0L) LE Table(n-1L)				;ascending flag
 IF sz(0) EQ 0 THEN BEGIN
	jl = -1L
	ju = n
 ENDIF ELSE BEGIN
	sz = sz(1:sz(0))
	jl = MAKE_ARRAY(DIMENSION=sz, VALUE=-1L)	;initial lower bounds
	ju = MAKE_ARRAY(DIMENSION=sz, VALUE=n)		;initial upper bounds
 ENDELSE
 WHILE (MAX(ju-jl) GT 1L) DO BEGIN			;binary search loop
	jm = (jl + ju)/2L				;midpoints
	IF asc THEN temp = Value GE Table(jm) ELSE temp = VALUE LE Table(jm)
	test = WHERE(temp, ntest)
	IF (ntest GT 0) THEN jl(test) = jm(test)
	test = WHERE(temp EQ 0, ntest)
	IF (ntest GT 0) THEN ju(test) = jm(test)
 ENDWHILE
 RETURN, jl
END
