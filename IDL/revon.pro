;--------------------------------------------------------------------
; RevOn Procedure (Changes to reversed cursor on VT100 Screen)
;--------------------------------------------------------------------
PRO RevOn
b=bytarr(4)
b(0)='1B'XB
b(1)=byte('[')
b(2)=byte('7')
b(3)=byte('m')
a=string(b)
PRINT,FORMAT='($,A4)',a
END
