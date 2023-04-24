using System;
using UnityEngine;

public class DynamicColor
{
	public enum E_BLINKING_MODE
    {
        BLINKING_NONE,
        BLINKING_BINARY,
        BLINKING_REGULAR,
        BLINKING_SMOOTH
    }

    // Pre-computed colors - colors "2" have 60% of opacity unlike the original (100%)
    public enum E_COLOR
    {
        RED,
        GREEN,
        BLUE,
        YELLOW,
        CYAN,
        MAGENTA,
        LIGHT_RED,
        LIGHT_GREEN,
        LIGHT_BLUE,
        LIGHT_CYAN,
        RED_2,
        GREEN_2,
        BLUE_2,
        YELLOW_2,
        CYAN_2,
        MAGENTA_2,
        LIGHT_RED_2,
        LIGHT_GREEN_2,
        LIGHT_BLUE_2,
        LIGHT_CYAN_2
    }

    private float _blinking_period;
    private E_BLINKING_MODE _blinking_type;
    private Color _dark_color;
    private Color _light_color;
    private float _time;

    public DynamicColor(Color lightColor, Color darkColor, E_BLINKING_MODE blinkingType, float blinkingPeriod)
    {
        _dark_color = darkColor;
        _light_color = lightColor;
        _blinking_type = blinkingType;
        _blinking_period = Mathf.Max(0.1f, blinkingPeriod);
        _time = Time.time;
    }

    public DynamicColor(E_COLOR color, E_BLINKING_MODE blinkingType, float blinkingPeriod)
    {
        _dark_color = darkColorsArray[(int)color];
        _light_color = lightColorsArray[(int)color];
        _blinking_type = blinkingType;
        _blinking_period = Mathf.Max(0.1f, blinkingPeriod);
        _time = Time.time;
    }

    public void changeBlinkingMode(E_BLINKING_MODE blinkingType, float blinkingPeriod)
    {
        _blinking_type = blinkingType;
        _blinking_period = Mathf.Max(0.1f, blinkingPeriod);
    }

    public void changeColor(Color lightColor, Color darkColor)
    {
        _light_color = lightColor;
        _dark_color = darkColor;
    }

    public void changeColor(E_COLOR color)
    {
        _light_color = lightColorsArray[(int)color];
        _dark_color = darkColorsArray[(int)color];
    }

    public Color getColor()
    {
        // case of fixed color
        if (_blinking_type == E_BLINKING_MODE.BLINKING_NONE) return _light_color;

        // Otherwise apply the blinking process
        float currentTime = Time.time;

        while(currentTime > _time + _blinking_period)
        {
            _time += _blinking_period;
        }

        float currentFractionOfPeriod = Mathf.Clamp01((currentTime - _time) / _blinking_period);

        float colorIntensity = calulateColorIntensity(currentFractionOfPeriod);

        return Color.Lerp(_dark_color, _light_color, colorIntensity);
    }

    public void restartBlinkingProcess()
    {
        _time = Time.time;
    }

    private float calulateColorIntensity(float theCurrentFractionOfPeriod)
    {
        float currentFractionOfPeriod = Mathf.Clamp01(theCurrentFractionOfPeriod);

        float colorIntensity = 0.0f;

        switch (_blinking_type)
        {
            case E_BLINKING_MODE.BLINKING_REGULAR:
                if (currentFractionOfPeriod < 0.5f)
                {
                    colorIntensity = 1.0f - 2.0f * currentFractionOfPeriod;
                }
                else
                {
                    colorIntensity = 2.0f * currentFractionOfPeriod - 1.0f;
                }
                break;

            case E_BLINKING_MODE.BLINKING_SMOOTH:
                float x;
                
                if(currentFractionOfPeriod > 0.5f)
                {
                    x = currentFractionOfPeriod - 0.5f;
                }
                else
                {
                    x = currentFractionOfPeriod + 0.5f;
                }
                colorIntensity = 4.0f * x * (1.0f - x);
                break;

            case E_BLINKING_MODE.BLINKING_BINARY:
                if(currentFractionOfPeriod < 0.5f)
                {
                    colorIntensity = 1.0f;
                }
                else
                {
                    colorIntensity = 0.0f;
                }
                break;

            case E_BLINKING_MODE.BLINKING_NONE:
            default:
                colorIntensity = 1.0f;
                break;
        }

        return colorIntensity;
    }

    private Color[] darkColorsArray = { new Color(0.400f, 0.000f, 0.000f),   // RED
                                        new Color(0.000f, 0.400f, 0.000f),   // GREEN
                                        new Color(0.000f, 0.000f, 0.400f),   // BLUE
                                        new Color(0.400f, 0.400f, 0.000f),   // YELLOW
                                        new Color(0.000f, 0.400f, 0.400f),   // CYAN
                                        new Color(0.400f, 0.000f, 0.400f),   // MAGENTA
                                        new Color(0.627f, 0.078f, 0.078f),   // LIGHT RED
                                        new Color(0.118f, 0.471f, 0.275f),   // LIGHT GREEN
                                        new Color(0.118f, 0.196f, 0.471f),   // LIGHT BLUE
                                        new Color(0.176f, 0.549f, 0.549f),   // LIGHT CYAN
                                        new Color(0.360f, 0.000f, 0.000f, 0.6f),   // RED 2
                                        new Color(0.000f, 0.400f, 0.000f, 0.6f),   // GREEN 2
                                        new Color(0.000f, 0.000f, 0.400f, 0.6f),   // BLUE 2
                                        new Color(0.400f, 0.400f, 0.000f, 0.6f),   // YELLOW 2
                                        new Color(0.000f, 0.400f, 0.400f, 0.6f),   // CYAN 2
                                        new Color(0.400f, 0.000f, 0.400f, 0.6f),   // MAGENTA 2
                                        new Color(0.627f, 0.078f, 0.078f, 0.6f),   // LIGHT RED 2
                                        new Color(0.118f, 0.471f, 0.275f, 0.6f),   // LIGHT GREEN 2
                                        new Color(0.118f, 0.196f, 0.471f, 0.6f),   // LIGHT BLUE 2
                                        new Color(0.176f, 0.549f, 0.549f, 0.6f) }; // LIGHT CYAN 2

private Color[] lightColorsArray = { new Color(1.000f, 0.000f, 0.000f),   // RED
                                         new Color(0.000f, 1.000f, 0.000f),   // GREEN
                                         new Color(0.000f, 0.000f, 1.000f),   // BLUE
                                         new Color(1.000f, 1.000f, 0.000f),   // YELLOW
                                         new Color(0.000f, 1.000f, 1.000f),   // CYAN
                                         new Color(1.000f, 0.000f, 1.000f),   // MAGENTA
                                         new Color(1.000f, 0.235f, 0.235f),   // LIGHT RED
                                         new Color(0.314f, 1.000f, 0.627f),   // LIGHT GREEN
                                         new Color(0.549f, 0.863f, 1.000f),   // LIGHT BLUE
                                         new Color(0.314f, 1.000f, 1.000f),   // LIGHT CYAN
                                         new Color(1.000f, 0.000f, 0.000f, 0.6f),   // RED 2
                                         new Color(0.000f, 1.000f, 0.000f, 0.6f),   // GREEN 2
                                         new Color(0.000f, 0.000f, 1.000f, 0.6f),   // BLUE 2
                                         new Color(1.000f, 1.000f, 0.000f, 0.6f),   // YELLOW 2
                                         new Color(0.000f, 1.000f, 1.000f, 0.6f),   // CYAN 2
                                         new Color(1.000f, 0.000f, 1.000f, 0.6f),   // MAGENTA 2
                                         new Color(1.000f, 0.235f, 0.235f, 0.6f),   // LIGHT RED 2
                                         new Color(0.314f, 1.000f, 0.627f, 0.6f),   // LIGHT GREEN 2
                                         new Color(0.549f, 0.863f, 1.000f, 0.6f),   // LIGHT BLUE 2
                                         new Color(0.314f, 1.000f, 1.000f, 0.6f) }; // LIGHT CYAN 2
}
