function [residual, g1, g2] = FA_static(y, x, params)
%
% Status : Computes static model for Dynare
%
% Inputs : 
%   y         [M_.endo_nbr by 1] double    vector of endogenous variables in declaration order
%   x         [M_.exo_nbr by 1] double     vector of exogenous variables in declaration order
%   params    [M_.param_nbr by 1] double   vector of parameter values in declaration order
%
% Outputs:
%   residual  [M_.endo_nbr by 1] double    vector of residuals of the static model equations 
%                                          in order of declaration of the equations
%   g1        [M_.endo_nbr by M_.endo_nbr] double    Jacobian matrix of the static model equations;
%                                                    columns: equations in order of declaration
%                                                    rows: variables in declaration order
%   g2        [M_.endo_nbr by (M_.endo_nbr)^2] double   Hessian matrix of the static model equations;
%                                                       columns: equations in order of declaration
%                                                       rows: variables in declaration order
%
%
% Warning : this file is generated automatically by Dynare
%           from model file (.mod)

residual = zeros( 39, 1);

%
% Model equations
%

T46 = exp(y(10))*exp(y(22))*(1-params(7))*exp(y(1))/exp(y(5));
T64 = exp(y(11))*params(1)*(1-params(6))*(exp(y(12))-exp(y(13)))+exp(y(17))*params(6)*params(1)*exp(y(11))*exp(y(21));
T122 = exp(y(22))*params(7)*exp(y(2))/exp(y(3))+exp(y(36))*(exp(y(9))-exp(y(32)));
T132 = exp(y(35))*(exp(y(3))*exp(y(36))*exp(y(25)))^params(7);
T152 = exp(y(22))*params(7)*exp(y(2))/exp(y(25));
T156 = exp(y(3))*exp(y(36))*params(31)*exp(y(25))^params(5);
T184 = exp(y(38));
T190 = exp(y(27))*params(11)*T184^((-params(12))*params(10));
T191 = T184^params(10);
T192 = T190*T191;
T201 = (1-params(11)*T184^(params(12)*(1-params(11)))*T184^(params(11)-1))/(1-params(11));
T221 = exp(y(28))*T191*exp(y(11))*params(1)*params(11)*T184^(params(12)*(-params(10)));
T231 = T184^(params(12)*(1-params(10)));
T233 = exp(y(29))*exp(y(11))*params(1)*params(11)*T184^(params(10)-1)*T231;
T241 = T184*exp(y(28))*params(10)/(params(10)-1)/exp(y(29));
T254 = exp(y(30))^params(15);
T258 = 1/params(1)*T184^params(13);
T259 = exp(y(26))/(params(10)/(params(10)-1));
T264 = (T258*T259^params(14))^(1-params(15));
T297 = exp(y(1))*exp(y(22))*params(7)/(exp(y(3))*exp(y(36)));
T317 = (-(exp(y(1))*exp(y(22))*(1-params(7))/exp(y(5))));
T320 = (-(exp(y(22))*params(7)*exp(y(2))/exp(y(3))/exp(y(9))));
T332 = (-(exp(y(5))^(1-params(7))*exp(y(35))*exp(y(3))*exp(y(36))*exp(y(25))*getPowerDeriv(exp(y(3))*exp(y(36))*exp(y(25)),params(7),1)));
T341 = (-((-(exp(y(3))*exp(y(36))*exp(y(1))*exp(y(22))*params(7)))/(exp(y(3))*exp(y(36))*exp(y(3))*exp(y(36)))));
T365 = (exp(y(7))-exp(y(7))*params(3))*getPowerDeriv(exp(y(7))-exp(y(7))*params(3),(-params(2)),1);
T436 = getPowerDeriv(T258*T259^params(14),1-params(15),1);
lhs =exp(y(10));
rhs =(exp(y(7))-exp(y(7))*params(3))^(-params(2))-(exp(y(7))-exp(y(7))*params(3))^(-params(2))*params(3)*params(1);
residual(1)= lhs-rhs;
lhs =params(1)*exp(y(13))*exp(y(11));
rhs =1;
residual(2)= lhs-rhs;
lhs =exp(y(11));
rhs =1;
residual(3)= lhs-rhs;
lhs =params(30)*exp(y(5))^params(4);
rhs =T46;
residual(4)= lhs-rhs;
lhs =exp(y(17));
rhs =T64;
residual(5)= lhs-rhs;
lhs =exp(y(18));
rhs =1-params(6)+exp(y(18))*params(6)*params(1)*exp(y(11))*exp(y(20));
residual(6)= lhs-rhs;
lhs =exp(y(19));
rhs =exp(y(18))/(params(29)-exp(y(17)));
residual(7)= lhs-rhs;
lhs =exp(y(20));
rhs =exp(y(13))+(exp(y(12))-exp(y(13)))*exp(y(19));
residual(8)= lhs-rhs;
lhs =exp(y(21));
rhs =exp(y(20));
residual(9)= lhs-rhs;
lhs =exp(y(9))*exp(y(3));
rhs =exp(y(19))*exp(y(14));
residual(10)= lhs-rhs;
lhs =exp(y(14));
rhs =exp(y(15))+exp(y(16));
residual(11)= lhs-rhs;
lhs =exp(y(15));
rhs =exp(y(14))*params(6)*exp(y(20))*exp((-x(4)));
residual(12)= lhs-rhs;
lhs =exp(y(16));
rhs =exp(y(3))*exp(y(9))*params(28)*exp(y(36));
residual(13)= lhs-rhs;
lhs =exp(y(12));
rhs =T122/exp(y(9));
residual(14)= lhs-rhs;
lhs =exp(y(2));
rhs =T132*exp(y(5))^(1-params(7));
residual(15)= lhs-rhs;
lhs =exp(y(9));
rhs =1;
residual(16)= lhs-rhs;
lhs =exp(y(32));
rhs =params(32)+params(31)/(1+params(5))*exp(y(25))^(1+params(5));
residual(17)= lhs-rhs;
lhs =T152;
rhs =T156;
residual(18)= lhs-rhs;
lhs =y(33);
rhs =exp(y(6))-exp(y(3))*exp(y(36))*exp(y(32));
residual(19)= lhs-rhs;
lhs =exp(y(3));
rhs =y(33)+exp(y(3))*exp(y(36));
residual(20)= lhs-rhs;
lhs =exp(y(8));
rhs =params(33)*exp(y(37));
residual(21)= lhs-rhs;
lhs =exp(y(1));
rhs =exp(y(6))+exp(y(7))+exp(y(8));
residual(22)= lhs-rhs;
lhs =exp(y(2));
rhs =exp(y(1))*exp(y(27));
residual(23)= lhs-rhs;
lhs =exp(y(27));
rhs =T192+(1-params(11))*T201^((-params(10))/(1-params(11)));
residual(24)= lhs-rhs;
lhs =exp(y(26));
rhs =1/exp(y(22));
residual(25)= lhs-rhs;
lhs =exp(y(28));
rhs =exp(y(22))*exp(y(1))+T221;
residual(26)= lhs-rhs;
lhs =exp(y(29));
rhs =exp(y(1))+T233;
residual(27)= lhs-rhs;
lhs =exp(y(39));
rhs =T241;
residual(28)= lhs-rhs;
lhs =T184^(1-params(10));
rhs =params(11)*T231+(1-params(11))*exp(y(39))^(1-params(10));
residual(29)= lhs-rhs;
lhs =exp(y(30));
rhs =exp(y(13))*T184;
residual(30)= lhs-rhs;
lhs =exp(y(30));
rhs =T254*T264*exp(x(5));
residual(31)= lhs-rhs;
lhs =y(35);
rhs =y(35)*params(18)-x(1);
residual(32)= lhs-rhs;
lhs =y(36);
rhs =y(36)*params(16)-x(2);
residual(33)= lhs-rhs;
lhs =y(37);
rhs =y(37)*params(20)-x(3);
residual(34)= lhs-rhs;
lhs =exp(y(4));
rhs =exp(y(3))*exp(y(36));
residual(35)= lhs-rhs;
lhs =exp(y(23));
rhs =exp(y(1))*exp(y(22))*(1-params(7))/exp(y(5));
residual(36)= lhs-rhs;
lhs =exp(y(24));
rhs =T297;
residual(37)= lhs-rhs;
lhs =y(34);
rhs =log(exp(y(7))-exp(y(7))*params(3))-params(30)*exp(y(5))^(1+params(4))/(1+params(4))+params(1)*y(34);
residual(38)= lhs-rhs;
lhs =exp(y(31));
rhs =exp(y(12))/exp(y(13));
residual(39)= lhs-rhs;
if ~isreal(residual)
  residual = real(residual)+imag(residual).^2;
end
if nargout >= 2,
  g1 = zeros(39, 39);

  %
  % Jacobian matrix
  %

  g1(1,7)=(-(T365-params(3)*params(1)*T365));
  g1(1,10)=exp(y(10));
  g1(2,11)=params(1)*exp(y(13))*exp(y(11));
  g1(2,13)=params(1)*exp(y(13))*exp(y(11));
  g1(3,11)=exp(y(11));
  g1(4,1)=(-T46);
  g1(4,5)=params(30)*exp(y(5))*getPowerDeriv(exp(y(5)),params(4),1)-(-(exp(y(5))*exp(y(10))*exp(y(22))*(1-params(7))*exp(y(1))))/(exp(y(5))*exp(y(5)));
  g1(4,10)=(-T46);
  g1(4,22)=(-T46);
  g1(5,11)=(-T64);
  g1(5,12)=(-(exp(y(11))*params(1)*(1-params(6))*exp(y(12))));
  g1(5,13)=(-(exp(y(11))*params(1)*(1-params(6))*(-exp(y(13)))));
  g1(5,17)=exp(y(17))-exp(y(17))*params(6)*params(1)*exp(y(11))*exp(y(21));
  g1(5,21)=(-(exp(y(17))*params(6)*params(1)*exp(y(11))*exp(y(21))));
  g1(6,11)=(-(exp(y(18))*params(6)*params(1)*exp(y(11))*exp(y(20))));
  g1(6,18)=exp(y(18))-exp(y(18))*params(6)*params(1)*exp(y(11))*exp(y(20));
  g1(6,20)=(-(exp(y(18))*params(6)*params(1)*exp(y(11))*exp(y(20))));
  g1(7,17)=(-((-(exp(y(18))*(-exp(y(17)))))/((params(29)-exp(y(17)))*(params(29)-exp(y(17))))));
  g1(7,18)=(-(exp(y(18))/(params(29)-exp(y(17)))));
  g1(7,19)=exp(y(19));
  g1(8,12)=(-(exp(y(12))*exp(y(19))));
  g1(8,13)=(-(exp(y(13))+exp(y(19))*(-exp(y(13)))));
  g1(8,19)=(-((exp(y(12))-exp(y(13)))*exp(y(19))));
  g1(8,20)=exp(y(20));
  g1(9,20)=(-exp(y(20)));
  g1(9,21)=exp(y(21));
  g1(10,3)=exp(y(9))*exp(y(3));
  g1(10,9)=exp(y(9))*exp(y(3));
  g1(10,14)=(-(exp(y(19))*exp(y(14))));
  g1(10,19)=(-(exp(y(19))*exp(y(14))));
  g1(11,14)=exp(y(14));
  g1(11,15)=(-exp(y(15)));
  g1(11,16)=(-exp(y(16)));
  g1(12,14)=(-(exp(y(14))*params(6)*exp(y(20))*exp((-x(4)))));
  g1(12,15)=exp(y(15));
  g1(12,20)=(-(exp(y(14))*params(6)*exp(y(20))*exp((-x(4)))));
  g1(13,3)=(-(exp(y(3))*exp(y(9))*params(28)*exp(y(36))));
  g1(13,9)=(-(exp(y(3))*exp(y(9))*params(28)*exp(y(36))));
  g1(13,16)=exp(y(16));
  g1(13,36)=(-(exp(y(3))*exp(y(9))*params(28)*exp(y(36))));
  g1(14,2)=T320;
  g1(14,3)=(-((-(exp(y(3))*exp(y(22))*params(7)*exp(y(2))))/(exp(y(3))*exp(y(3)))/exp(y(9))));
  g1(14,9)=(-((exp(y(9))*exp(y(9))*exp(y(36))-exp(y(9))*T122)/(exp(y(9))*exp(y(9)))));
  g1(14,12)=exp(y(12));
  g1(14,22)=T320;
  g1(14,32)=(-(exp(y(36))*(-exp(y(32)))/exp(y(9))));
  g1(14,36)=(-(exp(y(36))*(exp(y(9))-exp(y(32)))/exp(y(9))));
  g1(15,2)=exp(y(2));
  g1(15,3)=T332;
  g1(15,5)=(-(T132*exp(y(5))*getPowerDeriv(exp(y(5)),1-params(7),1)));
  g1(15,25)=T332;
  g1(15,35)=(-(T132*exp(y(5))^(1-params(7))));
  g1(15,36)=T332;
  g1(16,9)=exp(y(9));
  g1(17,25)=(-(params(31)/(1+params(5))*exp(y(25))*getPowerDeriv(exp(y(25)),1+params(5),1)));
  g1(17,32)=exp(y(32));
  g1(18,2)=T152;
  g1(18,3)=(-T156);
  g1(18,22)=T152;
  g1(18,25)=(-(exp(y(22))*params(7)*exp(y(2))*exp(y(25))))/(exp(y(25))*exp(y(25)))-exp(y(3))*exp(y(36))*params(31)*exp(y(25))*getPowerDeriv(exp(y(25)),params(5),1);
  g1(18,36)=(-T156);
  g1(19,3)=exp(y(3))*exp(y(36))*exp(y(32));
  g1(19,6)=(-exp(y(6)));
  g1(19,32)=exp(y(3))*exp(y(36))*exp(y(32));
  g1(19,33)=1;
  g1(19,36)=exp(y(3))*exp(y(36))*exp(y(32));
  g1(20,3)=exp(y(3))-exp(y(3))*exp(y(36));
  g1(20,33)=(-1);
  g1(20,36)=(-(exp(y(3))*exp(y(36))));
  g1(21,8)=exp(y(8));
  g1(21,37)=(-(params(33)*exp(y(37))));
  g1(22,1)=exp(y(1));
  g1(22,6)=(-exp(y(6)));
  g1(22,7)=(-exp(y(7)));
  g1(22,8)=(-exp(y(8)));
  g1(23,1)=(-(exp(y(1))*exp(y(27))));
  g1(23,2)=exp(y(2));
  g1(23,27)=(-(exp(y(1))*exp(y(27))));
  g1(24,27)=exp(y(27))-T192;
  g1(24,38)=(-(T191*exp(y(27))*params(11)*T184*getPowerDeriv(T184,(-params(12))*params(10),1)+T190*T184*getPowerDeriv(T184,params(10),1)+(1-params(11))*(-(T184^(params(11)-1)*params(11)*T184*getPowerDeriv(T184,params(12)*(1-params(11)),1)+params(11)*T184^(params(12)*(1-params(11)))*T184*getPowerDeriv(T184,params(11)-1,1)))/(1-params(11))*getPowerDeriv(T201,(-params(10))/(1-params(11)),1)));
  g1(25,22)=(-((-exp(y(22)))/(exp(y(22))*exp(y(22)))));
  g1(25,26)=exp(y(26));
  g1(26,1)=(-(exp(y(22))*exp(y(1))));
  g1(26,11)=(-T221);
  g1(26,22)=(-(exp(y(22))*exp(y(1))));
  g1(26,28)=exp(y(28))-T221;
  g1(26,38)=(-(exp(y(28))*(T184^(params(12)*(-params(10)))*exp(y(11))*params(1)*params(11)*T184*getPowerDeriv(T184,params(10),1)+T191*exp(y(11))*params(1)*params(11)*T184*getPowerDeriv(T184,params(12)*(-params(10)),1))));
  g1(27,1)=(-exp(y(1)));
  g1(27,11)=(-T233);
  g1(27,29)=exp(y(29))-T233;
  g1(27,38)=(-(exp(y(29))*(T231*exp(y(11))*params(1)*params(11)*T184*getPowerDeriv(T184,params(10)-1,1)+exp(y(11))*params(1)*params(11)*T184^(params(10)-1)*T184*getPowerDeriv(T184,params(12)*(1-params(10)),1))));
  g1(28,28)=(-T241);
  g1(28,29)=(-(T184*(-(exp(y(29))*exp(y(28))*params(10)/(params(10)-1)))/(exp(y(29))*exp(y(29)))));
  g1(28,38)=(-T241);
  g1(28,39)=exp(y(39));
  g1(29,38)=T184*getPowerDeriv(T184,1-params(10),1)-params(11)*T184*getPowerDeriv(T184,params(12)*(1-params(10)),1);
  g1(29,39)=(-((1-params(11))*exp(y(39))*getPowerDeriv(exp(y(39)),1-params(10),1)));
  g1(30,13)=(-(exp(y(13))*T184));
  g1(30,30)=exp(y(30));
  g1(30,38)=(-(exp(y(13))*T184));
  g1(31,26)=(-(exp(x(5))*T254*T258*T259*getPowerDeriv(T259,params(14),1)*T436));
  g1(31,30)=exp(y(30))-exp(x(5))*T264*exp(y(30))*getPowerDeriv(exp(y(30)),params(15),1);
  g1(31,38)=(-(exp(x(5))*T254*T436*T259^params(14)*1/params(1)*T184*getPowerDeriv(T184,params(13),1)));
  g1(32,35)=1-params(18);
  g1(33,36)=1-params(16);
  g1(34,37)=1-params(20);
  g1(35,3)=(-(exp(y(3))*exp(y(36))));
  g1(35,4)=exp(y(4));
  g1(35,36)=(-(exp(y(3))*exp(y(36))));
  g1(36,1)=T317;
  g1(36,5)=(-((-(exp(y(5))*exp(y(1))*exp(y(22))*(1-params(7))))/(exp(y(5))*exp(y(5)))));
  g1(36,22)=T317;
  g1(36,23)=exp(y(23));
  g1(37,1)=(-T297);
  g1(37,3)=T341;
  g1(37,22)=(-T297);
  g1(37,24)=exp(y(24));
  g1(37,36)=T341;
  g1(38,5)=params(30)*exp(y(5))*getPowerDeriv(exp(y(5)),1+params(4),1)/(1+params(4));
  g1(38,7)=(-1);
  g1(38,34)=1-params(1);
  g1(39,12)=(-(exp(y(12))/exp(y(13))));
  g1(39,13)=(-((-(exp(y(13))*exp(y(12))))/(exp(y(13))*exp(y(13)))));
  g1(39,31)=exp(y(31));
  if ~isreal(g1)
    g1 = real(g1)+2*imag(g1);
  end
end
if nargout >= 3,
  %
  % Hessian matrix
  %

  g2 = sparse([],[],[],39,1521);
end
end
