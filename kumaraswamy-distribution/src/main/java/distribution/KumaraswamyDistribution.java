package distribution;

import distribution.exceptions.NanException;
import distribution.exceptions.WrongIntervalException;
import distribution.handlers.Gamma;
import lombok.AllArgsConstructor;
import lombok.Getter;
import lombok.Setter;

@AllArgsConstructor
@Getter
@Setter
public class KumaraswamyDistribution {
    private Interval interval;


    /**
     * Probability Density Function, PDF
     * Функція густини ймовірності
     * функція густини ймовірності Кумарасвамі, розподіл без урахування інфляції
     * <p>
     * f(x;a,b) = a * b * x^(a-1) * (1-x^a)^(b-1)
     * <p>
     * Обмеження:
     * x є (0;1)
     * a,b > 0
     */

    public double pdf(double x) throws WrongIntervalException {
        if (x > 0 && x < 1 && interval.getA() > 0 && interval.getB() > 0)
            return interval.getA() * interval.getB() * Math.pow(x, interval.getA() - 1)
                    * Math.pow(1 - Math.pow(x, interval.getA()), interval.getB() - 1);
        else throw new WrongIntervalException();
    }

    /**
     * Cumulative distribution function, CDF
     * Кумулятивна функція розподілу
     * <p>
     * F(x;a,b) = ... = 1-(1-x^a)^b
     * <p>
     * Обмеження: x є (0;1)
     * TODO: УЗНАТЬ ЗА ДОПОЛНИТЕЛЬНЫЕ ОГРАНИЧЕНИЯ
     */

    public double cdf(double x) throws WrongIntervalException {
        if (x > 0 && x < 1)
            return 1 - Math.pow(1 - Math.pow(x, interval.getA()), interval.getB());
        else throw new WrongIntervalException();
    }


    /**
     * Percent point function, PPF
     * Inverse cumulative distribution function (Quantile function)
     * Зворотна кумулятивна функція розподілу (функція Квантилю)
     * <p>
     * F(y;a,b)^-1 = (1-(1-y)^1/b)^1/a
     * <p>
     * Обмеження: x є (0;1)
     * TODO: УЗНАТЬ ЗА ДОПОЛНИТЕЛЬНЫЕ ОГРАНИЧЕНИЯ
     */

    public double ppf(double x) throws WrongIntervalException {
        if (x > 0 && x < 1)
            return Math.pow(1 - Math.pow(1 - x, 1 / interval.getB()), 1 / interval.getA());
        else throw new WrongIntervalException();
    }

    /**
     * Survival function, SF
     * Функція виживання
     * <p>
     * Функції виживання найчастіше використовуються в надійності та суміжних областях.
     * Функція виживання – це ймовірність того, що змінна приймає значення більше за x.
     * <p>
     * S(x) = 1 - [1-(1-x^a)^b]
     * <p>
     * Обмеження: x є (0;1)
     * TODO: УЗНАТЬ ЗА ДОПОЛНИТЕЛЬНЫЕ ОГРАНИЧЕНИЯ
     */

    public double survivalFunction(double x) throws WrongIntervalException {
        if (x > 0 && x < 1)
            return 1 - cdf(x);
        else throw new WrongIntervalException();
    }

    /**
     * Hazard function, HF
     * Функція небезбечності
     * <p>
     * Функція небезпеки є відношенням функції щільності ймовірності (PDF)
     * до функції виживання (SF).
     * Діаграми небезпеки найчастіше використовуються в додатках надійності
     * <p>
     * Обмеження: x є (0;1)
     * <p>
     * TODO: УЗНАТЬ ЗА ДОПОЛНИТЕЛЬНЫЕ ОГРАНИЧЕНИЯ
     * TODO: УЗНАТЬ ЗА cumulative hazard function (СHF)
     * https://www.itl.nist.gov/div898/handbook/eda/section3/eda362.htm
     */

    public double hazardFunction(double x) throws WrongIntervalException {
        if (x > 0 && x < 1)
            return pdf(x) / survivalFunction(x);
        else throw new WrongIntervalException();
    }

    /**
     * medium, md
     * <p>
     * md = (1-2^(-1/b))1/a
     *  TODO: УЗНАТЬ ЗА ДОПОЛНИТЕЛЬНЫЕ ОГРАНИЧЕНИЯ
     */

    public double medium() {

        return Math.pow(1 - Math.pow(2, -1 / interval.getB()), 1 / interval.getA());
    }

    /**
     * mode
     * <p>
     * mode = (a-1/a*b-1)^1/a
     * <p>
     * Зауважте, що
     * - якщо a>1 і b>1 розподіл унімодальний
     * - якщо a<1 і b<1 розподіл уніантодальний
     * - якщо a>1 і b<=1 розподіл зростає
     * - якщо a<=1 і b>1 розподіл є спадним
     * - якщо a=b=1 розподіл постійний
     * Ця властивість визначена лише в перших двох випадках, повертаючи
     * mode та antimode відповідно.
     */

    public double mode() throws NanException {

        if ((interval.getA() > 1 && interval.getB() > 1) || (interval.getA() < 1 && interval.getB() < 1)) {
            return Math.pow(interval.getA() - 1 / interval.getA() * interval.getB() - 1, 1 / interval.getA());
        } else throw new NanException();
    }

    /**
     * Moments
     *
     * */
    public double moment(int n) {
        Gamma gamma = new Gamma();
        //Raw moment of order n
        // X=bΓ(1+1/a)Γ(b) / Γ(1+1a+b)*/
        return (interval.getB() * gamma.getGamma(1 + n / interval.getA()) *
                gamma.getGamma(interval.getB())) / gamma.getGamma(1 + n / interval.getA() + interval.getB());

    }

    /**
     * dispersion
     * Дисперсія в квадраті
     *
     * Дисперсія може бути обчислена з вихідних моментів
     * dispersion^2 = m2 - m1^2
     * */
    public double dispersionRaisedToTheSecondPower(){
        return moment(2) - Math.pow(moment(1),2);
    }

    /**
     * dispersion
     * Дисперсія
     *
     * Дисперсія може бути обчислена з вихідних моментів
     * dispersion = sqrt(m2 - m1^2)
     * */
    public double dispersion(){
        return Math.sqrt(dispersionRaisedToTheSecondPower());
    }
}

