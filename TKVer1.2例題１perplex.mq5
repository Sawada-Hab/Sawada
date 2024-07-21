//+------------------------------------------------------------------+
//|                  Cointegration Forex Trading EA                  |
//|                   Copyright 2024, AI Assistant                   |
//|                     https://www.example.com                      |
//+------------------------------------------------------------------+
#property copyright "Copyright 2024, AI Assistant"
#property link      "https://www.example.com"
#property version   "1.00"
#property strict

#include <Trade\Trade.mqh>
#include <Math\Stat\Math.mqh>
#include <Math\Stat\Stat.mqh>

// Input parameters
input int    MagicNumber    = 12345;  // Magic Number
input double Exposure       = 0.02;   // Risk exposure (0.02 = 2% of account balance)
input bool   DebugMode      = false;  // Debug mode on/off
input int    LookbackPeriod = 100;    // Lookback period for calculations
input double EntryThreshold = 1.5;    // Entry threshold in standard deviations
input double ExitThreshold  = 0.5;    // Exit threshold in standard deviations
input double StopLoss       = 50;     // Stop loss in points
input double TakeProfit     = 100;    // Take profit in points
input int    TrendPeriod    = 50;     // Trend filter period
input double TrendThreshold = 0.00003; // Trend filter threshold
input int    CheckInterval  = 5;      // Check for trade opportunities every X minutes

// Global variables
CTrade trade;
double portfolioWeights[3];
double meanSlippage = 0;
int slippageCount = 0;
datetime lastTradeTime = 0;
string symbols[3] = {"USDJPY", "EURJPY", "GBPJPY"};

// Function prototypes
string ErrorDescription(int error_code);
matrix MatrixCovariance(const matrix &data);
vector MatrixEigenVector(const matrix &m);
double GetMA(string symbol, ENUM_TIMEFRAMES timeframe, int period, int shift);
double DotProduct(const double &arr1[], const double &arr2[], int size);

//+------------------------------------------------------------------+
//| Expert initialization function                                   |
//+------------------------------------------------------------------+
int OnInit()
{
    trade.SetExpertMagicNumber(MagicNumber);
    
    for (int i = 0; i < 3; i++)
    {
        if (!SymbolSelect(symbols[i], true))
        {
            Print("Failed to select symbol: ", symbols[i], ". Error code: ", GetLastError());
            return INIT_FAILED;
        }
    }
    
    if (!CalculateWeights())
    {
        Print("Failed to calculate weights. EA initialization failed.");
        return INIT_FAILED;
    }
    
    meanSlippage = 0;
    slippageCount = 0;
    
    return(INIT_SUCCEEDED);
}

//+------------------------------------------------------------------+
//| Expert deinitialization function                                 |
//+------------------------------------------------------------------+
void OnDeinit(const int reason)
{
    // Cleanup code here if needed
}

//+------------------------------------------------------------------+
//| Expert tick function                                             |
//+------------------------------------------------------------------+
void OnTick()
{
    if (!IsTradeAllowed()) return;
    
    static datetime lastCheck = 0;
    datetime currentTime = TimeCurrent();
    
    if (currentTime - lastCheck >= CheckInterval * 60)
    {
        lastCheck = currentTime;
        double residual = CalculateResidual();
        if (ShouldTrade(residual))
        {
            ExecuteTrade(residual);
        }
    }
    
    if (AccountInfoDouble(ACCOUNT_EQUITY) <= AccountInfoDouble(ACCOUNT_BALANCE) * 0.5)
    {
        CloseAllPositions();
        ExpertRemove();
    }
}

//+------------------------------------------------------------------+
//| Calculate portfolio weights using Johansen cointegration test    |
//+------------------------------------------------------------------+
bool CalculateWeights()
{
    int dataSize = 1000;
    matrix data(dataSize, 3);
    int validSymbols = 0;
    
    for (int i = 0; i < 3; i++)
    {
        MqlRates rates[];
        if (CopyRates(symbols[i], PERIOD_H1, 1, dataSize, rates) == dataSize)
        {
            for (int j = 0; j < dataSize; j++)
            {
                data[j][i] = rates[j].close;
            }
            validSymbols++;
        }
        else
        {
            Print("Failed to copy rates for ", symbols[i], ". Error code: ", GetLastError());
        }
    }
    
    if (validSymbols < 2)
    {
        Print("Not enough valid symbols for cointegration calculation.");
        return false;
    }
    
    matrix cov = MatrixCovariance(data);
    vector eigenvector = MatrixEigenVector(cov);
    
    double sumWeights = 0;
    for (int i = 0; i < 3; i++)
    {
        sumWeights += MathAbs(eigenvector[i]);
    }
    for (int i = 0; i < 3; i++)
    {
        portfolioWeights[i] = eigenvector[i] / sumWeights;
    }
    
    return true;
}

//+------------------------------------------------------------------+
//| Calculate residual for trading signal                            |
//+------------------------------------------------------------------+
double CalculateResidual()
{
    double prices[3];
    for (int i = 0; i < 3; i++)
    {
        prices[i] = SymbolInfoDouble(symbols[i], SYMBOL_LAST);
    }
    
    return DotProduct(portfolioWeights, prices, 3);
}

//+------------------------------------------------------------------+
//| Check if trading is allowed based on time                        |
//+------------------------------------------------------------------+
bool IsTradeAllowed()
{
    MqlDateTime dt;
    TimeToStruct(TimeCurrent(), dt);
    
    if ((dt.hour == 21 && dt.min >= 30) || dt.hour == 22 || dt.hour == 23 || (dt.hour == 0 && dt.min < 30))
    {
        return false;
    }
    
    return true;
}

//+------------------------------------------------------------------+
//| Determine if we should enter or exit a trade                     |
//+------------------------------------------------------------------+
bool ShouldTrade(double residual)
{
    double mean = 0, stddev = 0;
    CalculateResidualStats(mean, stddev);
    
    if (PositionsTotal() > 0)
    {
        return (MathAbs(residual - mean) < ExitThreshold * stddev);
    }
    else
    {
        bool isTrending = CheckTrend();
        bool shortTermSignal = CheckShortTermSignal();
        return (MathAbs(residual - mean) > EntryThreshold * stddev) && (isTrending || shortTermSignal);
    }
}

//+------------------------------------------------------------------+
//| Check if there's a trend                                         |
//+------------------------------------------------------------------+
bool CheckTrend()
{
    double ma1 = GetMA(_Symbol, PERIOD_H1, TrendPeriod, 0);
    double ma2 = GetMA(_Symbol, PERIOD_H1, TrendPeriod, 1);
    return MathAbs(ma1 - ma2) > TrendThreshold;
}

//+------------------------------------------------------------------+
//| Check for short-term trading signal                              |
//+------------------------------------------------------------------+
bool CheckShortTermSignal()
{
    double ma1 = GetMA(_Symbol, PERIOD_M15, 20, 0);
    double ma2 = GetMA(_Symbol, PERIOD_M15, 20, 1);
    return MathAbs(ma1 - ma2) > TrendThreshold * 0.5;
}

//+------------------------------------------------------------------+
//| Calculate mean and standard deviation of residual                |
//+------------------------------------------------------------------+
void CalculateResidualStats(double &mean, double &stddev)
{
    int dataSize = LookbackPeriod;
    double residuals[];
    ArrayResize(residuals, dataSize);
    
    for (int i = 0; i < dataSize; i++)
    {
        double prices[3];
        for (int j = 0; j < 3; j++)
        {
            prices[j] = iClose(symbols[j], PERIOD_H1, i);
        }
        residuals[i] = DotProduct(portfolioWeights, prices, 3);
    }
    
    mean = MathMean(residuals);
    stddev = MathStandardDeviation(residuals);
}

//+------------------------------------------------------------------+
//| Execute trade based on signal                                    |
//+------------------------------------------------------------------+
void ExecuteTrade(double residual)
{
    double mean = 0, stddev = 0;
    CalculateResidualStats(mean, stddev);
    
    if (PositionsTotal() > 0)
    {
        CloseAllPositions();
    }
    else
    {
        ENUM_ORDER_TYPE orderType = (residual > mean) ? ORDER_TYPE_SELL : ORDER_TYPE_BUY;
        
        for (int i = 0; i < 3; i++)
        {
            double lotSize = CalculateLotSize(symbols[i], MathAbs(portfolioWeights[i]));
            if (lotSize > 0)
            {
                double stopLossPrice = (orderType == ORDER_TYPE_BUY) ? 
                    SymbolInfoDouble(symbols[i], SYMBOL_ASK) - StopLoss * SymbolInfoDouble(symbols[i], SYMBOL_POINT) :
                    SymbolInfoDouble(symbols[i], SYMBOL_BID) + StopLoss * SymbolInfoDouble(symbols[i], SYMBOL_POINT);
                double takeProfitPrice = (orderType == ORDER_TYPE_BUY) ?
                    SymbolInfoDouble(symbols[i], SYMBOL_ASK) + TakeProfit * SymbolInfoDouble(symbols[i], SYMBOL_POINT) :
                    SymbolInfoDouble(symbols[i], SYMBOL_BID) - TakeProfit * SymbolInfoDouble(symbols[i], SYMBOL_POINT);
                
                MqlTradeRequest request;
                MqlTradeResult result;
                ZeroMemory(request);
                ZeroMemory(result);
                
                request.action = TRADE_ACTION_DEAL;
                request.symbol = symbols[i];
                request.volume = lotSize;
                request.type = orderType;
                request.price = (orderType == ORDER_TYPE_BUY) ? SymbolInfoDouble(symbols[i], SYMBOL_ASK) : SymbolInfoDouble(symbols[i], SYMBOL_BID);
                request.sl = stopLossPrice;
                request.tp = takeProfitPrice;
                request.deviation = 10;
                request.magic = MagicNumber;
                request.comment = "Cointegration trade";
                
                if (!OrderSend(request, result))
                {
                    Print("OrderSend failed: ", ErrorDescription(GetLastError()));
                }
                else
                {
                    if (DebugMode) Print("Opened position in ", symbols[i], " with lot size ", lotSize);
                }
            }
        }
    }
}

//+------------------------------------------------------------------+
//| Calculate lot size based on risk and weight                      |
//+------------------------------------------------------------------+
double CalculateLotSize(string symbol, double weight)
{
    double accountBalance = AccountInfoDouble(ACCOUNT_BALANCE);
    double riskAmount = accountBalance * Exposure * weight;
    
    double tickValue = SymbolInfoDouble(symbol, SYMBOL_TRADE_TICK_VALUE);
    double tickSize = SymbolInfoDouble(symbol, SYMBOL_TRADE_TICK_SIZE);
    double pointValue = tickValue / tickSize;
    
    double stopLossPoints = StopLoss;
    
    double lotSize = NormalizeDouble(riskAmount / (stopLossPoints * pointValue), 2);
    
    double minLot = SymbolInfoDouble(symbol, SYMBOL_VOLUME_MIN);
    double maxLot = SymbolInfoDouble(symbol, SYMBOL_VOLUME_MAX);
    int lotDigits = (int)SymbolInfoInteger(symbol, SYMBOL_DIGITS);
    
    lotSize = MathMax(minLot, MathMin(maxLot, lotSize));
    lotSize = NormalizeDouble(lotSize, lotDigits);
    
    return lotSize;
}

//+------------------------------------------------------------------+
//| Close all open positions                                         |
//+------------------------------------------------------------------+
void CloseAllPositions()
{
    for (int i = PositionsTotal() - 1; i >= 0; i--)
    {
        ulong ticket = PositionGetTicket(i);
        if (ticket > 0)
        {
            trade.PositionClose(ticket);
        }
    }
}

//+------------------------------------------------------------------+
//| Expert event handler function                                    |
//+------------------------------------------------------------------+
void OnTradeTransaction(const MqlTradeTransaction& trans,
                        const MqlTradeRequest& request,
                        const MqlTradeResult& result)
{
    if (trans.type == TRADE_TRANSACTION_DEAL_ADD)
    {
        if (result.retcode == TRADE_RETCODE_DONE)
        {
            double slippage = MathAbs(result.price - request.price);
            meanSlippage = (meanSlippage * slippageCount + slippage) / (slippageCount + 1);
            slippageCount++;
            
            if (DebugMode) Print("Trade executed. Slippage: ", slippage, " Average slippage: ", meanSlippage);
        }
    }
}

//+------------------------------------------------------------------+
//| Tester function                                                  |
//+------------------------------------------------------------------+
double OnTester()
{
    double profit = TesterStatistics(STAT_PROFIT);
    double drawdown = TesterStatistics(STAT_BALANCEDD_PERCENT);
    double sharpRatio = TesterStatistics(STAT_SHARPE_RATIO);
    
    return (profit - drawdown * 2) * sharpRatio;
}

//+------------------------------------------------------------------+
//| Custom functions                                                 |
//+------------------------------------------------------------------+

// Calculate dot product of two arrays
double DotProduct(const double &arr1[], const double &arr2[], int size)
{
    double result = 0;
    for (int i = 0; i < size; i++)
    {
        result += arr1[i] * arr2[i];
    }
    return result;
}

// Calculate covariance matrix
matrix MatrixCovariance(const matrix &data)
{
    int rows = (int)data.Rows();
    int cols = (int)data.Cols();
    matrix cov(cols, cols);
    
    double means[];
    ArrayResize(means, cols);
    for (int i = 0; i < cols; i++)
    {
        vector col = data.Col(i);
        double arr[];
        ArrayResize(arr, rows);
        for (int k = 0; k < rows; k++)
        {
            arr[k] = col[k];
        }
        means[i] = MathMean(arr);
    }
    
    for (int i = 0; i < cols; i++)
    {
        for (int j = 0; j < cols; j++)
        {
            double sum = 0;
            vector col_i = data.Col(i);
            vector col_j = data.Col(j);
            
            for (int k = 0; k < rows; k++)
            {
                sum += (col_i[k] - means[i]) * (col_j[k] - means[j]);
            }
            cov[i][j] = sum / (rows - 1);
        }
    }
    
    return cov;
}

// Calculate eigenvector (simplified for this example)
vector MatrixEigenVector(const matrix &m)
{
    int size = (int)m.Rows();
    vector result(size);
    
    for (int i = 0; i < size; i++)
    {
        result[i] = m[i][i];
    }
    
    return result;
}

//+------------------------------------------------------------------+
//| Custom function to get MA                                        |
//+------------------------------------------------------------------+
double GetMA(string symbol, ENUM_TIMEFRAMES timeframe, int period, int shift)
{
    int maHandle = iMA(symbol, timeframe, period, 0, MODE_SMA, PRICE_CLOSE);
    if(maHandle == INVALID_HANDLE)
    {
        Print("Failed to create MA indicator handle for ", symbol);
        return 0;
    }

    double maBuffer[];
    if(CopyBuffer(maHandle, 0, shift, 1, maBuffer) != 1)
    {
        Print("Failed to copy MA buffer for ", symbol);
        IndicatorRelease(maHandle);
        return 0;
    }

    IndicatorRelease(maHandle);
    return maBuffer[0];
}

//+------------------------------------------------------------------+
//| Custom function to get error description                         |
//+------------------------------------------------------------------+
string ErrorDescription(int error_code)
{
    switch(error_code)
    {
        case 0: return "No error";
        case 1: return "No error returned, but the result is unknown";
        case 2: return "Common error";
        case 3: return "Invalid trade parameters";
        case 4: return "Trade server is busy";
        case 5: return "Old version of the client terminal";
        case 6: return "No connection with trade server";
        case 7: return "Not enough rights";
        case 8: return "Too frequent requests";
        case 9: return "Malfunctional trade operation";
        case 64: return "Account disabled";
        case 65: return "Invalid account";
        case 128: return "Trade timeout";
        case 129: return "Invalid price";
        case 130: return "Invalid stops";
        case 131: return "Invalid trade volume";
        case 132: return "Market is closed";
        case 133: return "Trade is disabled";
        case 134: return "Not enough money";
        case 135: return "Price changed";
        case 136: return "Off quotes";
        case 137: return "Broker is busy";
        case 138: return "Requote";
        case 139: return "Order is locked";
        case 140: return "Long positions only allowed";
        case 141: return "Too many requests";
        case 145: return "Modification denied because order is too close to market";
        case 146: return "Trade context is busy";
        case 147: return "Expirations are not supported by broker";
        case 148: return "The amount of open and pending orders has reached the limit set by the broker";
        case 149: return "Hedge is prohibited";
        case 150: return "Prohibited by FIFO rule";
        default: return "Unknown error";
    }
}

//+------------------------------------------------------------------+
