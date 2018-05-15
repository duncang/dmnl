function [DataOut] = SendData(DataIn)


Counter = DataIn(1);

Counter = Counter + 1;
DataInData = DataIn(Counter);

DataOut(2) = DataInData;

DataOut(1) = Counter;
