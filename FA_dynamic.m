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

residual = zeros(42, 1);
T40 = (1-params(1))*y(54)^((1-params(35))/params(36))+params(1)*y(55)^(1/params(36));
T87 = exp(y(24))*exp(y(36))*(1-params(7))*exp(y(15))/exp(y(19));
T107 = exp(y(58))*params(1)*(1-params(6))*(exp(y(59))-exp(y(27)))+params(6)*params(1)*exp(y(58))*exp(y(63))*exp(y(60));
T184 = exp(y(36))*params(7)*exp(y(16))/exp(y(1))+exp(y(50))*(exp(y(23))-exp(y(46)));
T196 = exp(y(49))*(exp(y(1))*exp(y(50))*exp(y(39)))^params(7);
T208 = (y(47)+params(34))/(params(34)+y(10))-1;
T210 = params(9)/2*T208^2;
T221 = params(1)*exp(y(58))*params(9)*((params(34)+y(66))/(y(47)+params(34))-1);
T222 = ((params(34)+y(66))/(y(47)+params(34)))^2;
T223 = T221*T222;
T235 = exp(y(36))*params(7)*exp(y(16))/exp(y(39));
T239 = exp(y(1))*exp(y(50))*params(31)*exp(y(39))^params(5);
T277 = params(11)*exp(y(8))*exp(y(14))^((-params(12))*params(10));
T281 = T277*exp(y(52))^params(10);
T290 = (1-params(11)*exp(y(14))^(params(12)*(1-params(11)))*exp(y(52))^(params(11)-1))/(1-params(11));
T315 = exp(y(58))*params(1)*params(11)*exp(y(68))^params(10)*exp(y(52))^(params(12)*(-params(10)))*exp(y(64));
T322 = exp(y(58))*params(1)*params(11)*exp(y(68))^(params(10)-1);
T329 = T322*exp(y(52))^(params(12)*(1-params(10)))*exp(y(65));
T337 = exp(y(52))*exp(y(42))*params(10)/(params(10)-1)/exp(y(43));
T353 = exp(y(9))^params(15);
T357 = 1/params(1)*exp(y(52))^params(13);
T358 = exp(y(40))/(params(10)/(params(10)-1));
T363 = (T357*T358^params(14))^(1-params(15));
T399 = exp(y(15))*exp(y(36))*params(7)/(exp(y(50))*exp(y(1)));
T419 = (-(exp(y(15))*exp(y(36))*(1-params(7))/exp(y(19))));
T422 = (-(exp(y(36))*params(7)*exp(y(16))/exp(y(1))/exp(y(3))));
T434 = (-(exp(y(19))^(1-params(7))*exp(y(49))*exp(y(1))*exp(y(50))*exp(y(39))*getPowerDeriv(exp(y(1))*exp(y(50))*exp(y(39)),params(7),1)));
T442 = (-((-(exp(y(50))*exp(y(1))*exp(y(15))*exp(y(36))*params(7)))/(exp(y(50))*exp(y(1))*exp(y(50))*exp(y(1)))));
T571 = getPowerDeriv(T357*T358^params(14),1-params(15),1);
T598 = params(9)/2*(-(y(47)+params(34)))/((params(34)+y(10))*(params(34)+y(10)))*2*T208;
T610 = params(9)/2*2*T208*1/(params(34)+y(10));
T655 = getPowerDeriv(T290,(-params(10))/(1-params(11)),1);
T714 = getPowerDeriv(T40,params(36)/(1-params(35)),1);
lhs =y(54);
rhs =log(y(21)-params(3)*y(2))-params(30)/(1+params(4))*y(19)^(1+params(4));
residual(1)= lhs-rhs;
lhs =y(55);
rhs =y(69)^(1-params(35));
residual(2)= lhs-rhs;
lhs =y(56);
rhs =T40^(params(36)/(1-params(35)));
residual(3)= lhs-rhs;
lhs =exp(y(24));
rhs =(exp(y(21))-params(3)*exp(y(2)))^(-params(2))-params(3)*params(1)*(exp(y(57))-params(3)*exp(y(21)))^(-params(2));
residual(4)= lhs-rhs;
lhs =params(1)*exp(y(27))*exp(y(58));
rhs =1;
residual(5)= lhs-rhs;
lhs =exp(y(25));
rhs =exp(y(24))/exp(y(4));
residual(6)= lhs-rhs;
lhs =params(30)*exp(y(19))^params(4);
rhs =T87;
residual(7)= lhs-rhs;
lhs =exp(y(31));
rhs =T107;
residual(8)= lhs-rhs;
lhs =exp(y(32));
rhs =1-params(6)+params(6)*params(1)*exp(y(58))*exp(y(62))*exp(y(61));
residual(9)= lhs-rhs;
lhs =exp(y(33));
rhs =exp(y(32))/(params(29)-exp(y(31)));
residual(10)= lhs-rhs;
lhs =exp(y(34));
rhs =exp(y(5))+(exp(y(26))-exp(y(5)))*exp(y(7));
residual(11)= lhs-rhs;
lhs =exp(y(35));
rhs =exp(y(34))*exp(y(33))/exp(y(7));
residual(12)= lhs-rhs;
lhs =exp(y(23))*exp(y(17));
rhs =exp(y(33))*exp(y(28));
residual(13)= lhs-rhs;
lhs =exp(y(28));
rhs =exp(y(29))+exp(y(30));
residual(14)= lhs-rhs;
lhs =exp(y(29));
rhs =params(6)*exp(y(34))*exp(y(6))*exp((-x(it_, 4)));
residual(15)= lhs-rhs;
lhs =exp(y(30));
rhs =exp(y(23))*params(28)*exp(y(50))*exp(y(1));
residual(16)= lhs-rhs;
lhs =exp(y(26));
rhs =T184/exp(y(3));
residual(17)= lhs-rhs;
lhs =exp(y(16));
rhs =T196*exp(y(19))^(1-params(7));
residual(18)= lhs-rhs;
lhs =exp(y(23));
rhs =1+T210+(y(47)+params(34))*params(9)*T208/(params(34)+y(10))-T223;
residual(19)= lhs-rhs;
lhs =exp(y(46));
rhs =params(32)+params(31)/(1+params(5))*exp(y(39))^(1+params(5));
residual(20)= lhs-rhs;
lhs =T235;
rhs =T239;
residual(21)= lhs-rhs;
lhs =y(47);
rhs =exp(y(20))-exp(y(1))*exp(y(50))*exp(y(46));
residual(22)= lhs-rhs;
lhs =exp(y(17));
rhs =y(47)+exp(y(50))*exp(y(1));
residual(23)= lhs-rhs;
lhs =exp(y(22));
rhs =params(33)*exp(y(51));
residual(24)= lhs-rhs;
lhs =exp(y(15));
rhs =exp(y(20))+exp(y(21))+exp(y(22))+(y(47)+params(34))*T210;
residual(25)= lhs-rhs;
lhs =exp(y(16));
rhs =exp(y(15))*exp(y(41));
residual(26)= lhs-rhs;
lhs =exp(y(41));
rhs =T281+(1-params(11))*T290^((-params(10))/(1-params(11)));
residual(27)= lhs-rhs;
lhs =exp(y(40));
rhs =1/exp(y(36));
residual(28)= lhs-rhs;
lhs =exp(y(42));
rhs =exp(y(36))*exp(y(15))+T315;
residual(29)= lhs-rhs;
lhs =exp(y(43));
rhs =exp(y(15))+T329;
residual(30)= lhs-rhs;
lhs =exp(y(53));
rhs =T337;
residual(31)= lhs-rhs;
lhs =exp(y(52))^(1-params(10));
rhs =params(11)*exp(y(14))^(params(12)*(1-params(10)))+(1-params(11))*exp(y(53))^(1-params(10));
residual(32)= lhs-rhs;
lhs =exp(y(44));
rhs =exp(y(27))*exp(y(68));
residual(33)= lhs-rhs;
lhs =exp(y(44));
rhs =T353*T363*exp(x(it_, 5));
residual(34)= lhs-rhs;
lhs =y(49);
rhs =params(18)*y(11)-x(it_, 1);
residual(35)= lhs-rhs;
lhs =y(50);
rhs =params(16)*y(12)-x(it_, 2);
residual(36)= lhs-rhs;
lhs =y(51);
rhs =params(20)*y(13)-x(it_, 3);
residual(37)= lhs-rhs;
lhs =exp(y(18));
rhs =exp(y(50))*exp(y(1));
residual(38)= lhs-rhs;
lhs =exp(y(37));
rhs =exp(y(15))*exp(y(36))*(1-params(7))/exp(y(19));
residual(39)= lhs-rhs;
lhs =exp(y(38));
rhs =T399;
residual(40)= lhs-rhs;
lhs =y(48);
rhs =log(exp(y(21))-params(3)*exp(y(2)))-params(30)*exp(y(19))^(1+params(4))/(1+params(4))+params(1)*y(67);
residual(41)= lhs-rhs;
lhs =exp(y(45));
rhs =exp(y(59))/exp(y(27));
residual(42)= lhs-rhs;
if nargout >= 2,
  g1 = zeros(42, 74);

  %
  % Jacobian matrix
  %

  g1(1,19)=params(30)/(1+params(4))*getPowerDeriv(y(19),1+params(4),1);
  g1(1,2)=(-((-params(3))/(y(21)-params(3)*y(2))));
  g1(1,21)=(-(1/(y(21)-params(3)*y(2))));
  g1(1,54)=1;
  g1(2,55)=1;
  g1(2,69)=(-(getPowerDeriv(y(69),1-params(35),1)));
  g1(3,54)=(-((1-params(1))*getPowerDeriv(y(54),(1-params(35))/params(36),1)*T714));
  g1(3,55)=(-(T714*params(1)*getPowerDeriv(y(55),1/params(36),1)));
  g1(3,56)=1;
  g1(4,2)=(-((-(params(3)*exp(y(2))))*getPowerDeriv(exp(y(21))-params(3)*exp(y(2)),(-params(2)),1)));
  g1(4,21)=(-(exp(y(21))*getPowerDeriv(exp(y(21))-params(3)*exp(y(2)),(-params(2)),1)-params(3)*params(1)*(-(params(3)*exp(y(21))))*getPowerDeriv(exp(y(57))-params(3)*exp(y(21)),(-params(2)),1)));
  g1(4,57)=params(3)*params(1)*exp(y(57))*getPowerDeriv(exp(y(57))-params(3)*exp(y(21)),(-params(2)),1);
  g1(4,24)=exp(y(24));
  g1(5,58)=params(1)*exp(y(27))*exp(y(58));
  g1(5,27)=params(1)*exp(y(27))*exp(y(58));
  g1(6,4)=(-((-(exp(y(24))*exp(y(4))))/(exp(y(4))*exp(y(4)))));
  g1(6,24)=(-(exp(y(24))/exp(y(4))));
  g1(6,25)=exp(y(25));
  g1(7,15)=(-T87);
  g1(7,19)=params(30)*exp(y(19))*getPowerDeriv(exp(y(19)),params(4),1)-(-(exp(y(19))*exp(y(24))*exp(y(36))*(1-params(7))*exp(y(15))))/(exp(y(19))*exp(y(19)));
  g1(7,24)=(-T87);
  g1(7,36)=(-T87);
  g1(8,58)=(-T107);
  g1(8,59)=(-(exp(y(58))*params(1)*(1-params(6))*exp(y(59))));
  g1(8,27)=(-(exp(y(58))*params(1)*(1-params(6))*(-exp(y(27)))));
  g1(8,31)=exp(y(31));
  g1(8,60)=(-(params(6)*params(1)*exp(y(58))*exp(y(63))*exp(y(60))));
  g1(8,63)=(-(params(6)*params(1)*exp(y(58))*exp(y(63))*exp(y(60))));
  g1(9,58)=(-(params(6)*params(1)*exp(y(58))*exp(y(62))*exp(y(61))));
  g1(9,32)=exp(y(32));
  g1(9,61)=(-(params(6)*params(1)*exp(y(58))*exp(y(62))*exp(y(61))));
  g1(9,62)=(-(params(6)*params(1)*exp(y(58))*exp(y(62))*exp(y(61))));
  g1(10,31)=(-((-(exp(y(32))*(-exp(y(31)))))/((params(29)-exp(y(31)))*(params(29)-exp(y(31))))));
  g1(10,32)=(-(exp(y(32))/(params(29)-exp(y(31)))));
  g1(10,33)=exp(y(33));
  g1(11,26)=(-(exp(y(26))*exp(y(7))));
  g1(11,5)=(-(exp(y(5))+exp(y(7))*(-exp(y(5)))));
  g1(11,7)=(-((exp(y(26))-exp(y(5)))*exp(y(7))));
  g1(11,34)=exp(y(34));
  g1(12,7)=(-(exp(y(34))*(-(exp(y(33))*exp(y(7))))/(exp(y(7))*exp(y(7)))));
  g1(12,33)=(-(exp(y(34))*exp(y(33))/exp(y(7))));
  g1(12,34)=(-(exp(y(34))*exp(y(33))/exp(y(7))));
  g1(12,35)=exp(y(35));
  g1(13,17)=exp(y(23))*exp(y(17));
  g1(13,23)=exp(y(23))*exp(y(17));
  g1(13,28)=(-(exp(y(33))*exp(y(28))));
  g1(13,33)=(-(exp(y(33))*exp(y(28))));
  g1(14,28)=exp(y(28));
  g1(14,29)=(-exp(y(29)));
  g1(14,30)=(-exp(y(30)));
  g1(15,6)=(-(params(6)*exp(y(34))*exp(y(6))*exp((-x(it_, 4)))));
  g1(15,29)=exp(y(29));
  g1(15,34)=(-(params(6)*exp(y(34))*exp(y(6))*exp((-x(it_, 4)))));
  g1(15,73)=(-(params(6)*exp(y(34))*exp(y(6))*(-exp((-x(it_, 4))))));
  g1(16,1)=(-(exp(y(23))*params(28)*exp(y(50))*exp(y(1))));
  g1(16,23)=(-(exp(y(23))*params(28)*exp(y(50))*exp(y(1))));
  g1(16,30)=exp(y(30));
  g1(16,50)=(-(exp(y(23))*params(28)*exp(y(50))*exp(y(1))));
  g1(17,16)=T422;
  g1(17,1)=(-((-(exp(y(1))*exp(y(36))*params(7)*exp(y(16))))/(exp(y(1))*exp(y(1)))/exp(y(3))));
  g1(17,3)=(-((-(T184*exp(y(3))))/(exp(y(3))*exp(y(3)))));
  g1(17,23)=(-(exp(y(23))*exp(y(50))/exp(y(3))));
  g1(17,26)=exp(y(26));
  g1(17,36)=T422;
  g1(17,46)=(-(exp(y(50))*(-exp(y(46)))/exp(y(3))));
  g1(17,50)=(-(exp(y(50))*(exp(y(23))-exp(y(46)))/exp(y(3))));
  g1(18,16)=exp(y(16));
  g1(18,1)=T434;
  g1(18,19)=(-(T196*exp(y(19))*getPowerDeriv(exp(y(19)),1-params(7),1)));
  g1(18,39)=T434;
  g1(18,49)=(-(T196*exp(y(19))^(1-params(7))));
  g1(18,50)=T434;
  g1(19,23)=exp(y(23));
  g1(19,58)=T223;
  g1(19,10)=(-(T598+((params(34)+y(10))*(y(47)+params(34))*params(9)*(-(y(47)+params(34)))/((params(34)+y(10))*(params(34)+y(10)))-(y(47)+params(34))*params(9)*T208)/((params(34)+y(10))*(params(34)+y(10)))));
  g1(19,47)=(-(T610+(params(9)*T208+(y(47)+params(34))*params(9)*1/(params(34)+y(10)))/(params(34)+y(10))-(T222*params(1)*exp(y(58))*params(9)*(-(params(34)+y(66)))/((y(47)+params(34))*(y(47)+params(34)))+T221*(-(params(34)+y(66)))/((y(47)+params(34))*(y(47)+params(34)))*2*(params(34)+y(66))/(y(47)+params(34)))));
  g1(19,66)=T222*params(1)*exp(y(58))*params(9)*1/(y(47)+params(34))+T221*2*(params(34)+y(66))/(y(47)+params(34))*1/(y(47)+params(34));
  g1(20,39)=(-(params(31)/(1+params(5))*exp(y(39))*getPowerDeriv(exp(y(39)),1+params(5),1)));
  g1(20,46)=exp(y(46));
  g1(21,16)=T235;
  g1(21,1)=(-T239);
  g1(21,36)=T235;
  g1(21,39)=(-(exp(y(36))*params(7)*exp(y(16))*exp(y(39))))/(exp(y(39))*exp(y(39)))-exp(y(1))*exp(y(50))*params(31)*exp(y(39))*getPowerDeriv(exp(y(39)),params(5),1);
  g1(21,50)=(-T239);
  g1(22,1)=exp(y(1))*exp(y(50))*exp(y(46));
  g1(22,20)=(-exp(y(20)));
  g1(22,46)=exp(y(1))*exp(y(50))*exp(y(46));
  g1(22,47)=1;
  g1(22,50)=exp(y(1))*exp(y(50))*exp(y(46));
  g1(23,1)=(-(exp(y(50))*exp(y(1))));
  g1(23,17)=exp(y(17));
  g1(23,47)=(-1);
  g1(23,50)=(-(exp(y(50))*exp(y(1))));
  g1(24,22)=exp(y(22));
  g1(24,51)=(-(params(33)*exp(y(51))));
  g1(25,15)=exp(y(15));
  g1(25,20)=(-exp(y(20)));
  g1(25,21)=(-exp(y(21)));
  g1(25,22)=(-exp(y(22)));
  g1(25,10)=(-((y(47)+params(34))*T598));
  g1(25,47)=(-(T210+(y(47)+params(34))*T610));
  g1(26,15)=(-(exp(y(15))*exp(y(41))));
  g1(26,16)=exp(y(16));
  g1(26,41)=(-(exp(y(15))*exp(y(41))));
  g1(27,8)=(-T281);
  g1(27,41)=exp(y(41));
  g1(27,14)=(-(exp(y(52))^params(10)*params(11)*exp(y(8))*exp(y(14))*getPowerDeriv(exp(y(14)),(-params(12))*params(10),1)+(1-params(11))*(-(exp(y(52))^(params(11)-1)*params(11)*exp(y(14))*getPowerDeriv(exp(y(14)),params(12)*(1-params(11)),1)))/(1-params(11))*T655));
  g1(27,52)=(-(T277*exp(y(52))*getPowerDeriv(exp(y(52)),params(10),1)+(1-params(11))*T655*(-(params(11)*exp(y(14))^(params(12)*(1-params(11)))*exp(y(52))*getPowerDeriv(exp(y(52)),params(11)-1,1)))/(1-params(11))));
  g1(28,36)=(-((-exp(y(36)))/(exp(y(36))*exp(y(36)))));
  g1(28,40)=exp(y(40));
  g1(29,15)=(-(exp(y(36))*exp(y(15))));
  g1(29,58)=(-T315);
  g1(29,36)=(-(exp(y(36))*exp(y(15))));
  g1(29,42)=exp(y(42));
  g1(29,64)=(-T315);
  g1(29,52)=(-(exp(y(64))*exp(y(58))*params(1)*params(11)*exp(y(68))^params(10)*exp(y(52))*getPowerDeriv(exp(y(52)),params(12)*(-params(10)),1)));
  g1(29,68)=(-(exp(y(64))*exp(y(52))^(params(12)*(-params(10)))*exp(y(58))*params(1)*params(11)*exp(y(68))*getPowerDeriv(exp(y(68)),params(10),1)));
  g1(30,15)=(-exp(y(15)));
  g1(30,58)=(-T329);
  g1(30,43)=exp(y(43));
  g1(30,65)=(-T329);
  g1(30,52)=(-(exp(y(65))*T322*exp(y(52))*getPowerDeriv(exp(y(52)),params(12)*(1-params(10)),1)));
  g1(30,68)=(-(exp(y(65))*exp(y(52))^(params(12)*(1-params(10)))*exp(y(58))*params(1)*params(11)*exp(y(68))*getPowerDeriv(exp(y(68)),params(10)-1,1)));
  g1(31,42)=(-T337);
  g1(31,43)=(-(exp(y(52))*(-(exp(y(43))*exp(y(42))*params(10)/(params(10)-1)))/(exp(y(43))*exp(y(43)))));
  g1(31,52)=(-T337);
  g1(31,53)=exp(y(53));
  g1(32,14)=(-(params(11)*exp(y(14))*getPowerDeriv(exp(y(14)),params(12)*(1-params(10)),1)));
  g1(32,52)=exp(y(52))*getPowerDeriv(exp(y(52)),1-params(10),1);
  g1(32,53)=(-((1-params(11))*exp(y(53))*getPowerDeriv(exp(y(53)),1-params(10),1)));
  g1(33,27)=(-(exp(y(27))*exp(y(68))));
  g1(33,44)=exp(y(44));
  g1(33,68)=(-(exp(y(27))*exp(y(68))));
  g1(34,40)=(-(exp(x(it_, 5))*T353*T357*T358*getPowerDeriv(T358,params(14),1)*T571));
  g1(34,9)=(-(exp(x(it_, 5))*T363*exp(y(9))*getPowerDeriv(exp(y(9)),params(15),1)));
  g1(34,44)=exp(y(44));
  g1(34,52)=(-(exp(x(it_, 5))*T353*T571*T358^params(14)*1/params(1)*exp(y(52))*getPowerDeriv(exp(y(52)),params(13),1)));
  g1(34,74)=(-(T353*T363*exp(x(it_, 5))));
  g1(35,11)=(-params(18));
  g1(35,49)=1;
  g1(35,70)=1;
  g1(36,12)=(-params(16));
  g1(36,50)=1;
  g1(36,71)=1;
  g1(37,13)=(-params(20));
  g1(37,51)=1;
  g1(37,72)=1;
  g1(38,1)=(-(exp(y(50))*exp(y(1))));
  g1(38,18)=exp(y(18));
  g1(38,50)=(-(exp(y(50))*exp(y(1))));
  g1(39,15)=T419;
  g1(39,19)=(-((-(exp(y(19))*exp(y(15))*exp(y(36))*(1-params(7))))/(exp(y(19))*exp(y(19)))));
  g1(39,36)=T419;
  g1(39,37)=exp(y(37));
  g1(40,15)=(-T399);
  g1(40,1)=T442;
  g1(40,36)=(-T399);
  g1(40,38)=exp(y(38));
  g1(40,50)=T442;
  g1(41,19)=params(30)*exp(y(19))*getPowerDeriv(exp(y(19)),1+params(4),1)/(1+params(4));
  g1(41,2)=(-((-(params(3)*exp(y(2))))/(exp(y(21))-params(3)*exp(y(2)))));
  g1(41,21)=(-(exp(y(21))/(exp(y(21))-params(3)*exp(y(2)))));
  g1(41,48)=1;
  g1(41,67)=(-params(1));
  g1(42,59)=(-(exp(y(59))/exp(y(27))));
  g1(42,27)=(-((-(exp(y(27))*exp(y(59))))/(exp(y(27))*exp(y(27)))));
  g1(42,45)=exp(y(45));
end
if nargout >= 3,
  %
  % Hessian matrix
  %

  g2 = sparse([],[],[],42,5476);
end
if nargout >= 4,
  %
  % Third order derivatives
  %

  g3 = sparse([],[],[],42,405224);
end
end
