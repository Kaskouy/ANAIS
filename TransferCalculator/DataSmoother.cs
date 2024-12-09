using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using static SFS.Parsers.Ini.IniDataEnv;

class DataSmoother
{
    private readonly uint _nbMaxValues;

    private readonly List<double> _weights;

    private LinkedList<Double2> _values;

    // Creates a DataSmoother, designed to absorb small variations of a particular data.
    // This is done by memorizing several values over time, and by averaging all values.
    // All values can have a specific weight. For consistency, the most recent values have the most weight.
    // ---------------------------------------------------------------------------------
    // Parameters:
    // - nbMaxValues: the max number of values that will be memorized - must be greater than 1
    // - minWeight: the weight given to the oldest value
    // - maxWeight: the weight given to the newest value
    // - quadraticBehaviourFraction: indicates how "quadratic" the weight function is:
    //   * 0.0 : weight varies linearly
    //   * 1.0 : weight varies quadratically, with the minimum corresponding exactly to the minWeight value
    //   * values below 0.0 or higher than 1.0 will work, but would give weird results
    // ---------------------------------------------------------------------------------
    public DataSmoother(uint nbMaxValues, double minWeight, double maxWeight, double quadraticBehaviourFraction)
    {
        _nbMaxValues = nbMaxValues;
        _values = new LinkedList<Double2>();
        _values.Clear();

        _weights = new List<double>();
        InitWeightTable(minWeight, maxWeight, quadraticBehaviourFraction);
    }

    private void InitWeightTable(double minWeight, double maxWeight, double quadraticBehaviourFraction)
    {
        // Represents the length of the interval between 0 and (_nbMaxValues - 1), which are the indexes of the weight table.
        double deltaX = (double)(_nbMaxValues - 1);

        // weight will be easily expressed as weight = a * index^2 + b*index + c, where index = 0 for the first value, and index = _nbMaxValues - 1 for the last one.
        double a = quadraticBehaviourFraction * (maxWeight - minWeight) / (deltaX * deltaX);
        double b = (minWeight - maxWeight) / deltaX - a * (_nbMaxValues - 1);
        double c = maxWeight;

        for (uint index = 0; index < _nbMaxValues; index++)
        {
            _weights.Add(a * index * index + b * index + c);
        }
    }

    public void Add(Double2 value) 
    {
        // insert new element at the beginning of the list
        _values.AddFirst(value);

        // Remove the last value if we've reached the maximum number of values
        if (_values.Count > _nbMaxValues)
        {
            _values.RemoveLast();
        }
    }

    public void Clear()
    {
        _values.Clear();
    }

    public Double2 GetSmoothedData()
    {
        Double2 smoothedData = new Double2(0.0, 0.0);
        double weightSum = 0.0;
        int index = 0;

        foreach(Double2 value in _values)
        {
            smoothedData = smoothedData + value * _weights[index];
            weightSum += _weights[index];
            index++;
        }

        smoothedData /= weightSum;

        return smoothedData;
    }

}
