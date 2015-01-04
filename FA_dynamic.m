function [residual, g1, g2, g3] = FA_dynamic(y, x, params, steady_state, it_)
%
% Status : Computes dynamic model for Dynare
%
% Inputs :
%   y         [#dynamic variables by 1] double    vector of endogenous variables in the order stored
%                                                 in M_.lead_lag_incidence; see the Manual
%   x         [M_.exo_nbr by nperiods] double     matrix of exogenous variables (in declaration order)
%                                                 for all simulation periods
%   params    [M_.param_nbr by 1] double          vector of parameter values in declaration order
%   it_       scalar double                       time period for exogenous variables for which to evaluate the model
%
% Outputs:
%   residual  [M_.endo_nbr by 1] double    vector of residuals of the dynamic model equations in order of 
%                                          declaration of the equations
%   g1        [M_.endo_nbr by #dynamic variables] double    Jacobian matrix of the dynamic model equations;
%                                                           columns: equations in order of declaration
%                                                           rows: variables in order stored in M_.lead_lag_incidence
%   g2        [M_.endo_nbr by (#dynamic variables)^2] double   Hessian matrix of the dynamic model equations;
%                                                              columns: equations in order of declaration
%                                                              rows: variables in order stored in M_.lead_lag_incidence
%   g3        [M_.endo_nbr by (#dynamic variables)^3] double   Third order derivative matrix of the dynamic model equations;
%                                                              columns: equations in order of declaration
%                                                              rows: variables in order stored in M_.lead_lag_incidence
%
%
% Warning : this file is generated automatically by Dynare
%           from model file (.mod)

%
% Model equations
%

residual = zeros(39, 1);
T58 = exp(y(24))*exp(y(36))*(1-params(7))*exp(y(15))/exp(y(19));
T78 = exp(y(55))*params(1)*(1-params(6))*(exp(y(56))-exp(y(27)))+params(6)*params(1)*exp(y(55))*exp(y(60))*exp(y(57));
T155 = exp(y(36))*params(7)*exp(y(16))/exp(y(1))+exp(y(50))*(exp(y(23))-exp(y(46)));
T167 = exp(y(49))*(exp(y(1))*exp(y(50))*exp(y(39)))^params(7);
T179 = (y(47)+params(34))/(params(34)+y(10))-1;
T181 = params(9)/2*T179^2;
T192 = params(1)*exp(y(55))*params(9)*((params(34)+y(63))/(y(47)+params(34))-1);
T193 = ((params(34)+y(63))/(y(47)+params(34)))^2;
T194 = T192*T193;
T206 = exp(y(36))*params(7)*exp(y(16))/exp(y(39));
T210 = exp(y(1))*exp(y(50))*params(31)*exp(y(39))^params(5);
T248 = params(11)*exp(y(8))*exp(y(14))^((-params(12))*params(10));
T252 = T248*exp(y(52))^params(10);
T261 = (1-params(11)*exp(y(14))^(params(12)*(1-params(11)))*exp(y(52))^(params(11)-1))/(1-params(11));
T286 = exp(y(55))*params(1)*params(11)*exp(y(65))^params(10)*exp(y(52))^(params(12)*(-params(10)))*exp(y(61));
T293 = exp(y(55))*params(1)*params(11)*exp(y(65))^(params(10)-1);
T300 = T293*exp(y(52))^(params(12)*(1-params(10)))*exp(y(62));
T308 = exp(y(52))*exp(y(42))*params(10)/(params(10)-1)/exp(y(43));
T324 = exp(y(9))^params(15);
T328 = 1/params(1)*exp(y(52))^params(13);
T329 = exp(y(40))/(params(10)/(params(10)-1));
T334 = (T328*T329^params(14))^(1-params(15));
T370 = exp(y(15))*exp(y(36))*params(7)/(exp(y(50))*exp(y(1)));
T391 = (-(exp(y(15))*exp(y(36))*(1-params(7))/exp(y(19))));
T394 = (-(exp(y(36))*params(7)*exp(y(16))/exp(y(1))/exp(y(3))));
T406 = (-(exp(y(19))^(1-params(7))*exp(y(49))*exp(y(1))*exp(y(50))*exp(y(39))*getPowerDeriv(exp(y(1))*exp(y(50))*exp(y(39)),params(7),1)));
T414 = (-((-(exp(y(50))*exp(y(1))*exp(y(15))*exp(y(36))*params(7)))/(exp(y(50))*exp(y(1))*exp(y(50))*exp(y(1)))));
T535 = getPowerDeriv(T328*T329^params(14),1-params(15),1);
T562 = params(9)/2*(-(y(47)+params(34)))/((params(34)+y(10))*(params(34)+y(10)))*2*T179;
T574 = params(9)/2*2*T179*1/(params(34)+y(10));
T619 = getPowerDeriv(T261,(-params(10))/(1-params(11)),1);
lhs =exp(y(24));
rhs =(exp(y(21))-params(3)*exp(y(2)))^(-params(2))-params(3)*params(1)*(exp(y(54))-exp(y(21))*params(3))^(-params(2));
residual(1)= lhs-rhs;
lhs =params(1)*exp(y(27))*exp(y(55));
rhs =1;
residual(2)= lhs-rhs;
lhs =exp(y(25));
rhs =exp(y(24))/exp(y(4));
residual(3)= lhs-rhs;
lhs =params(30)*exp(y(19))^params(4);
rhs =T58;
residual(4)= lhs-rhs;
lhs =exp(y(31));
rhs =T78;
residual(5)= lhs-rhs;
lhs =exp(y(32));
rhs =1-params(6)+params(6)*params(1)*exp(y(55))*exp(y(59))*exp(y(58));
residual(6)= lhs-rhs;
lhs =exp(y(33));
rhs =exp(y(32))/(params(29)-exp(y(31)));
residual(7)= lhs-rhs;
lhs =exp(y(34));
rhs =exp(y(5))+(exp(y(26))-exp(y(5)))*exp(y(7));
residual(8)= lhs-rhs;
lhs =exp(y(35));
rhs =exp(y(34))*exp(y(33))/exp(y(7));
residual(9)= lhs-rhs;
lhs =exp(y(23))*exp(y(17));
rhs =exp(y(33))*exp(y(28));
residual(10)= lhs-rhs;
lhs =exp(y(28));
rhs =exp(y(29))+exp(y(30));
residual(11)= lhs-rhs;
lhs =exp(y(29));
rhs =params(6)*exp(y(34))*exp(y(6))*exp((-x(it_, 4)));
residual(12)= lhs-rhs;
lhs =exp(y(30));
rhs =exp(y(23))*params(28)*exp(y(50))*exp(y(1));
residual(13)= lhs-rhs;
lhs =exp(y(26));
rhs =T155/exp(y(3));
residual(14)= lhs-rhs;
lhs =exp(y(16));
rhs =T167*exp(y(19))^(1-params(7));
residual(15)= lhs-rhs;
lhs =exp(y(23));
rhs =1+T181+(y(47)+params(34))*params(9)*T179/(params(34)+y(10))-T194;
residual(16)= lhs-rhs;
lhs =exp(y(46));
rhs =params(32)+params(31)/(1+params(5))*exp(y(39))^(1+params(5));
residual(17)= lhs-rhs;
lhs =T206;
rhs =T210;
residual(18)= lhs-rhs;
lhs =y(47);
rhs =exp(y(20))-exp(y(1))*exp(y(50))*exp(y(46));
residual(19)= lhs-rhs;
lhs =exp(y(17));
rhs =y(47)+exp(y(50))*exp(y(1));
residual(20)= lhs-rhs;
lhs =exp(y(22));
rhs =params(33)*exp(y(51));
residual(21)= lhs-rhs;
lhs =exp(y(15));
rhs =exp(y(20))+exp(y(21))+exp(y(22))+(y(47)+params(34))*T181;
residual(22)= lhs-rhs;
lhs =exp(y(16));
rhs =exp(y(15))*exp(y(41));
residual(23)= lhs-rhs;
lhs =exp(y(41));
rhs =T252+(1-params(11))*T261^((-params(10))/(1-params(11)));
residual(24)= lhs-rhs;
lhs =exp(y(40));
rhs =1/exp(y(36));
residual(25)= lhs-rhs;
lhs =exp(y(42));
rhs =exp(y(36))*exp(y(15))+T286;
residual(26)= lhs-rhs;
lhs =exp(y(43));
rhs =exp(y(15))+T300;
residual(27)= lhs-rhs;
lhs =exp(y(53));
rhs =T308;
residual(28)= lhs-rhs;
lhs =exp(y(52))^(1-params(10));
rhs =params(11)*exp(y(14))^(params(12)*(1-params(10)))+(1-params(11))*exp(y(53))^(1-params(10));
residual(29)= lhs-rhs;
lhs =exp(y(44));
rhs =exp(y(27))*exp(y(65));
residual(30)= lhs-rhs;
lhs =exp(y(44));
rhs =T324*T334*exp(x(it_, 5));
residual(31)= lhs-rhs;
lhs =y(49);
rhs =params(18)*y(11)-x(it_, 1);
residual(32)= lhs-rhs;
lhs =y(50);
rhs =params(16)*y(12)-x(it_, 2);
residual(33)= lhs-rhs;
lhs =y(51);
rhs =params(20)*y(13)-x(it_, 3);
residual(34)= lhs-rhs;
lhs =exp(y(18));
rhs =exp(y(50))*exp(y(1));
residual(35)= lhs-rhs;
lhs =exp(y(37));
rhs =exp(y(15))*exp(y(36))*(1-params(7))/exp(y(19));
residual(36)= lhs-rhs;
lhs =exp(y(38));
rhs =T370;
residual(37)= lhs-rhs;
lhs =y(48);
rhs =log(exp(y(21))-params(3)*exp(y(2)))-params(30)*exp(y(19))^(1+params(4))/(1+params(4))+params(1)*y(64);
residual(38)= lhs-rhs;
lhs =exp(y(45));
rhs =exp(y(56))/exp(y(27));
residual(39)= lhs-rhs;
if nargout >= 2,
  g1 = zeros(39, 70);

  %
  % Jacobian matrix
  %

  g1(1,2)=(-((-(params(3)*exp(y(2))))*getPowerDeriv(exp(y(21))-params(3)*exp(y(2)),(-params(2)),1)));
  g1(1,21)=(-(exp(y(21))*getPowerDeriv(exp(y(21))-params(3)*exp(y(2)),(-params(2)),1)-params(3)*params(1)*(-(exp(y(21))*params(3)))*getPowerDeriv(exp(y(54))-exp(y(21))*params(3),(-params(2)),1)));
  g1(1,54)=params(3)*params(1)*exp(y(54))*getPowerDeriv(exp(y(54))-exp(y(21))*params(3),(-params(2)),1);
  g1(1,24)=exp(y(24));
  g1(2,55)=params(1)*exp(y(27))*exp(y(55));
  g1(2,27)=params(1)*exp(y(27))*exp(y(55));
  g1(3,4)=(-((-(exp(y(24))*exp(y(4))))/(exp(y(4))*exp(y(4)))));
  g1(3,24)=(-(exp(y(24))/exp(y(4))));
  g1(3,25)=exp(y(25));
  g1(4,15)=(-T58);
  g1(4,19)=params(30)*exp(y(19))*getPowerDeriv(exp(y(19)),params(4),1)-(-(exp(y(19))*exp(y(24))*exp(y(36))*(1-params(7))*exp(y(15))))/(exp(y(19))*exp(y(19)));
  g1(4,24)=(-T58);
  g1(4,36)=(-T58);
  g1(5,55)=(-T78);
  g1(5,56)=(-(exp(y(55))*params(1)*(1-params(6))*exp(y(56))));
  g1(5,27)=(-(exp(y(55))*params(1)*(1-params(6))*(-exp(y(27)))));
  g1(5,31)=exp(y(31));
  g1(5,57)=(-(params(6)*params(1)*exp(y(55))*exp(y(60))*exp(y(57))));
  g1(5,60)=(-(params(6)*params(1)*exp(y(55))*exp(y(60))*exp(y(57))));
  g1(6,55)=(-(params(6)*params(1)*exp(y(55))*exp(y(59))*exp(y(58))));
  g1(6,32)=exp(y(32));
  g1(6,58)=(-(params(6)*params(1)*exp(y(55))*exp(y(59))*exp(y(58))));
  g1(6,59)=(-(params(6)*params(1)*exp(y(55))*exp(y(59))*exp(y(58))));
  g1(7,31)=(-((-(exp(y(32))*(-exp(y(31)))))/((params(29)-exp(y(31)))*(params(29)-exp(y(31))))));
  g1(7,32)=(-(exp(y(32))/(params(29)-exp(y(31)))));
  g1(7,33)=exp(y(33));
  g1(8,26)=(-(exp(y(26))*exp(y(7))));
  g1(8,5)=(-(exp(y(5))+exp(y(7))*(-exp(y(5)))));
  g1(8,7)=(-((exp(y(26))-exp(y(5)))*exp(y(7))));
  g1(8,34)=exp(y(34));
  g1(9,7)=(-(exp(y(34))*(-(exp(y(33))*exp(y(7))))/(exp(y(7))*exp(y(7)))));
  g1(9,33)=(-(exp(y(34))*exp(y(33))/exp(y(7))));
  g1(9,34)=(-(exp(y(34))*exp(y(33))/exp(y(7))));
  g1(9,35)=exp(y(35));
  g1(10,17)=exp(y(23))*exp(y(17));
  g1(10,23)=exp(y(23))*exp(y(17));
  g1(10,28)=(-(exp(y(33))*exp(y(28))));
  g1(10,33)=(-(exp(y(33))*exp(y(28))));
  g1(11,28)=exp(y(28));
  g1(11,29)=(-exp(y(29)));
  g1(11,30)=(-exp(y(30)));
  g1(12,6)=(-(params(6)*exp(y(34))*exp(y(6))*exp((-x(it_, 4)))));
  g1(12,29)=exp(y(29));
  g1(12,34)=(-(params(6)*exp(y(34))*exp(y(6))*exp((-x(it_, 4)))));
  g1(12,69)=(-(params(6)*exp(y(34))*exp(y(6))*(-exp((-x(it_, 4))))));
  g1(13,1)=(-(exp(y(23))*params(28)*exp(y(50))*exp(y(1))));
  g1(13,23)=(-(exp(y(23))*params(28)*exp(y(50))*exp(y(1))));
  g1(13,30)=exp(y(30));
  g1(13,50)=(-(exp(y(23))*params(28)*exp(y(50))*exp(y(1))));
  g1(14,16)=T394;
  g1(14,1)=(-((-(exp(y(1))*exp(y(36))*params(7)*exp(y(16))))/(exp(y(1))*exp(y(1)))/exp(y(3))));
  g1(14,3)=(-((-(T155*exp(y(3))))/(exp(y(3))*exp(y(3)))));
  g1(14,23)=(-(exp(y(23))*exp(y(50))/exp(y(3))));
  g1(14,26)=exp(y(26));
  g1(14,36)=T394;
  g1(14,46)=(-(exp(y(50))*(-exp(y(46)))/exp(y(3))));
  g1(14,50)=(-(exp(y(50))*(exp(y(23))-exp(y(46)))/exp(y(3))));
  g1(15,16)=exp(y(16));
  g1(15,1)=T406;
  g1(15,19)=(-(T167*exp(y(19))*getPowerDeriv(exp(y(19)),1-params(7),1)));
  g1(15,39)=T406;
  g1(15,49)=(-(T167*exp(y(19))^(1-params(7))));
  g1(15,50)=T406;
  g1(16,23)=exp(y(23));
  g1(16,55)=T194;
  g1(16,10)=(-(T562+((params(34)+y(10))*(y(47)+params(34))*params(9)*(-(y(47)+params(34)))/((params(34)+y(10))*(params(34)+y(10)))-(y(47)+params(34))*params(9)*T179)/((params(34)+y(10))*(params(34)+y(10)))));
  g1(16,47)=(-(T574+(params(9)*T179+(y(47)+params(34))*params(9)*1/(params(34)+y(10)))/(params(34)+y(10))-(T193*params(1)*exp(y(55))*params(9)*(-(params(34)+y(63)))/((y(47)+params(34))*(y(47)+params(34)))+T192*(-(params(34)+y(63)))/((y(47)+params(34))*(y(47)+params(34)))*2*(params(34)+y(63))/(y(47)+params(34)))));
  g1(16,63)=T193*params(1)*exp(y(55))*params(9)*1/(y(47)+params(34))+T192*2*(params(34)+y(63))/(y(47)+params(34))*1/(y(47)+params(34));
  g1(17,39)=(-(params(31)/(1+params(5))*exp(y(39))*getPowerDeriv(exp(y(39)),1+params(5),1)));
  g1(17,46)=exp(y(46));
  g1(18,16)=T206;
  g1(18,1)=(-T210);
  g1(18,36)=T206;
  g1(18,39)=(-(exp(y(36))*params(7)*exp(y(16))*exp(y(39))))/(exp(y(39))*exp(y(39)))-exp(y(1))*exp(y(50))*params(31)*exp(y(39))*getPowerDeriv(exp(y(39)),params(5),1);
  g1(18,50)=(-T210);
  g1(19,1)=exp(y(1))*exp(y(50))*exp(y(46));
  g1(19,20)=(-exp(y(20)));
  g1(19,46)=exp(y(1))*exp(y(50))*exp(y(46));
  g1(19,47)=1;
  g1(19,50)=exp(y(1))*exp(y(50))*exp(y(46));
  g1(20,1)=(-(exp(y(50))*exp(y(1))));
  g1(20,17)=exp(y(17));
  g1(20,47)=(-1);
  g1(20,50)=(-(exp(y(50))*exp(y(1))));
  g1(21,22)=exp(y(22));
  g1(21,51)=(-(params(33)*exp(y(51))));
  g1(22,15)=exp(y(15));
  g1(22,20)=(-exp(y(20)));
  g1(22,21)=(-exp(y(21)));
  g1(22,22)=(-exp(y(22)));
  g1(22,10)=(-((y(47)+params(34))*T562));
  g1(22,47)=(-(T181+(y(47)+params(34))*T574));
  g1(23,15)=(-(exp(y(15))*exp(y(41))));
  g1(23,16)=exp(y(16));
  g1(23,41)=(-(exp(y(15))*exp(y(41))));
  g1(24,8)=(-T252);
  g1(24,41)=exp(y(41));
  g1(24,14)=(-(exp(y(52))^params(10)*params(11)*exp(y(8))*exp(y(14))*getPowerDeriv(exp(y(14)),(-params(12))*params(10),1)+(1-params(11))*(-(exp(y(52))^(params(11)-1)*params(11)*exp(y(14))*getPowerDeriv(exp(y(14)),params(12)*(1-params(11)),1)))/(1-params(11))*T619));
  g1(24,52)=(-(T248*exp(y(52))*getPowerDeriv(exp(y(52)),params(10),1)+(1-params(11))*T619*(-(params(11)*exp(y(14))^(params(12)*(1-params(11)))*exp(y(52))*getPowerDeriv(exp(y(52)),params(11)-1,1)))/(1-params(11))));
  g1(25,36)=(-((-exp(y(36)))/(exp(y(36))*exp(y(36)))));
  g1(25,40)=exp(y(40));
  g1(26,15)=(-(exp(y(36))*exp(y(15))));
  g1(26,55)=(-T286);
  g1(26,36)=(-(exp(y(36))*exp(y(15))));
  g1(26,42)=exp(y(42));
  g1(26,61)=(-T286);
  g1(26,52)=(-(exp(y(61))*exp(y(55))*params(1)*params(11)*exp(y(65))^params(10)*exp(y(52))*getPowerDeriv(exp(y(52)),params(12)*(-params(10)),1)));
  g1(26,65)=(-(exp(y(61))*exp(y(52))^(params(12)*(-params(10)))*exp(y(55))*params(1)*params(11)*exp(y(65))*getPowerDeriv(exp(y(65)),params(10),1)));
  g1(27,15)=(-exp(y(15)));
  g1(27,55)=(-T300);
  g1(27,43)=exp(y(43));
  g1(27,62)=(-T300);
  g1(27,52)=(-(exp(y(62))*T293*exp(y(52))*getPowerDeriv(exp(y(52)),params(12)*(1-params(10)),1)));
  g1(27,65)=(-(exp(y(62))*exp(y(52))^(params(12)*(1-params(10)))*exp(y(55))*params(1)*params(11)*exp(y(65))*getPowerDeriv(exp(y(65)),params(10)-1,1)));
  g1(28,42)=(-T308);
  g1(28,43)=(-(exp(y(52))*(-(exp(y(43))*exp(y(42))*params(10)/(params(10)-1)))/(exp(y(43))*exp(y(43)))));
  g1(28,52)=(-T308);
  g1(28,53)=exp(y(53));
  g1(29,14)=(-(params(11)*exp(y(14))*getPowerDeriv(exp(y(14)),params(12)*(1-params(10)),1)));
  g1(29,52)=exp(y(52))*getPowerDeriv(exp(y(52)),1-params(10),1);
  g1(29,53)=(-((1-params(11))*exp(y(53))*getPowerDeriv(exp(y(53)),1-params(10),1)));
  g1(30,27)=(-(exp(y(27))*exp(y(65))));
  g1(30,44)=exp(y(44));
  g1(30,65)=(-(exp(y(27))*exp(y(65))));
  g1(31,40)=(-(exp(x(it_, 5))*T324*T328*T329*getPowerDeriv(T329,params(14),1)*T535));
  g1(31,9)=(-(exp(x(it_, 5))*T334*exp(y(9))*getPowerDeriv(exp(y(9)),params(15),1)));
  g1(31,44)=exp(y(44));
  g1(31,52)=(-(exp(x(it_, 5))*T324*T535*T329^params(14)*1/params(1)*exp(y(52))*getPowerDeriv(exp(y(52)),params(13),1)));
  g1(31,70)=(-(T324*T334*exp(x(it_, 5))));
  g1(32,11)=(-params(18));
  g1(32,49)=1;
  g1(32,66)=1;
  g1(33,12)=(-params(16));
  g1(33,50)=1;
  g1(33,67)=1;
  g1(34,13)=(-params(20));
  g1(34,51)=1;
  g1(34,68)=1;
  g1(35,1)=(-(exp(y(50))*exp(y(1))));
  g1(35,18)=exp(y(18));
  g1(35,50)=(-(exp(y(50))*exp(y(1))));
  g1(36,15)=T391;
  g1(36,19)=(-((-(exp(y(19))*exp(y(15))*exp(y(36))*(1-params(7))))/(exp(y(19))*exp(y(19)))));
  g1(36,36)=T391;
  g1(36,37)=exp(y(37));
  g1(37,15)=(-T370);
  g1(37,1)=T414;
  g1(37,36)=(-T370);
  g1(37,38)=exp(y(38));
  g1(37,50)=T414;
  g1(38,19)=params(30)*exp(y(19))*getPowerDeriv(exp(y(19)),1+params(4),1)/(1+params(4));
  g1(38,2)=(-((-(params(3)*exp(y(2))))/(exp(y(21))-params(3)*exp(y(2)))));
  g1(38,21)=(-(exp(y(21))/(exp(y(21))-params(3)*exp(y(2)))));
  g1(38,48)=1;
  g1(38,64)=(-params(1));
  g1(39,56)=(-(exp(y(56))/exp(y(27))));
  g1(39,27)=(-((-(exp(y(27))*exp(y(56))))/(exp(y(27))*exp(y(27)))));
  g1(39,45)=exp(y(45));
end
if nargout >= 3,
  %
  % Hessian matrix
  %

  g2 = sparse([],[],[],39,4900);
end
if nargout >= 4,
  %
  % Third order derivatives
  %

  g3 = sparse([],[],[],39,343000);
end
end
