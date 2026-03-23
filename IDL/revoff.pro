;--------------------------------------------------------------------
; RevOff Procedure (Changes to regular cursor on VT100 Screen)
;--------------------------------------------------------------------
PRO RevOff
b=bytarr(4)
b(0)='1B'XB
b(1)=byte('[')
b(2)=byte('0')
b(3)=byte('m')
a=string(b)
PRINT,FORMAT='($,A4)',a
END
