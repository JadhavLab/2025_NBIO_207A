% poly3fit.m is a function that returns sse, the sum of squared errors which
% will be minimized with respect to parameters for a fitted 3rd order
% polynomial
function sse=poly4fit(params,Input,Actual_Output)
A0=params(1);
A1=params(2);
A2=params(3);
A3=params(4);

Fitted_Curve=A0+A1*Input+A2*Input.*Input+A3*Input.*Input.*Input;
Error_Vector=Fitted_Curve - Actual_Output;
% When curvefitting, a typical quantity to
% minimize is the sum of squares error
sse=sum(Error_Vector.^2);
% You could also write sse as
% sse=Error_Vector(:)'*Error_Vector(:);