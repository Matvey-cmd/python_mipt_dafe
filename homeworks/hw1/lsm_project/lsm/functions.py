"""
В этом модуле хранятся функции для применения МНК
"""

from typing import Optional
from numbers import Real  # раскомментируйте при необходимости

from lsm_project.event_logger.event_logger import EventLogger

from lsm_project.lsm.enumerations import MismatchStrategies
from lsm_project.lsm.models import (
    LSMDescription,
    LSMLines,
)

PRECISION = 3  # константа для точности вывода
event_logger = EventLogger()  # для логирования


def get_lsm_description(
        abscissa: list[float], ordinates: list[float],
        mismatch_strategy: MismatchStrategies = MismatchStrategies.FALL
) -> LSMDescription:
    """
    Функции для получения описания рассчитаной зависимости

    :param: abscissa - значения абсцисс
    :param: ordinates - значение ординат
    :param: mismatch_strategy - стратегия обработки несовпадения

    :return: структура типа LSMDescription
    """

    global event_logger

    _is_valid_measurments(abscissa)
    _is_valid_measurments(ordinates)
    if len(abscissa) != len(ordinates):
        _process_mismatch(abscissa, ordinates, mismatch_strategy)
    return _get_lsm_description(abscissa, ordinates)


def get_lsm_lines(
        abscissa: list[float], ordinates: list[float],
        lsm_description: Optional[LSMDescription] = None
) -> LSMLines:
    """
    Функция для расчета значений функций с помощью результатов МНК

    :param: abscissa - значения абсцисс
    :param: ordinates - значение ординат
    :param: lsm_description - описание МНК

    :return: структура типа LSMLines
    """

    a = lsm_description

    if a is None:
        a = get_lsm_description(abscissa, ordinates)
    if type(a) is not LSMDescription:
        raise TypeError
    y = [a.incline * x + a.shift for x in abscissa]
    above = [(
                         a.incline + a.incline_error) * x + a.shift + a.shift_error
             for x in abscissa]
    under = [(
                         a.incline - a.incline_error) * x + a.shift - a.shift_error
             for x in abscissa]

    return LSMLines(
        abscissa=abscissa,
        ordinates=ordinates,
        line_predicted=y,
        line_above=above,
        line_under=under
    )


def get_report(
        lsm_description: LSMDescription, path_to_save: str = ''
) -> str:
    """
    Функция для формирования отчета о результатах МНК

    :param: lsm_description - описание МНК
    :param: path_to_save - путь к файлу для сохранения отчета

    :return: строка - отчет определенного формата
    """
    global PRECISION

    report = '\n'.join([
        "="*40 + "LSM computing result" + "="*40 + "\n",
        "[INFO]: incline: " + f'{lsm_description.incline:.{PRECISION}f}' + ";",
        "[INFO]: shift: " + f'{lsm_description.shift:.{PRECISION}f}' + ";",
        "[INFO]: incline error: " + f'{lsm_description.incline_error:.{PRECISION}f}' + ";",
        "[INFO]: shift error: " + f'{lsm_description.shift_error:.{PRECISION}f}' + ";",
        "\n" + "="*100

    ])

    if path_to_save != '':
        with open(path_to_save, 'w') as f:
            f.write(report)
    return report


# служебная функция для валидации
def _is_valid_measurments(measurments: list[float]) -> bool:
    if len(measurments) <= 2:
        raise ValueError
    if not (all(isinstance(i, Real) for i in measurments)):
        raise ValueError
    if not (type(list(measurments)) is list):
        raise TypeError


# служебная функция для обработки несоответствия размеров
def _process_mismatch(
        abscissa: list[float], ordinates: list[float],
        mismatch_strategy: MismatchStrategies = MismatchStrategies.FALL
) -> tuple[list[float], list[float]]:
    global event_logger

    if mismatch_strategy == MismatchStrategies.FALL:
        raise RuntimeError
    elif mismatch_strategy == MismatchStrategies.CUT:
        if len(abscissa) < len(ordinates):
            while len(abscissa) < len(ordinates):
                ordinates.pop(-1)
        if len(abscissa) > len(ordinates):
            while len(abscissa) > len(ordinates):
                abscissa.pop(-1)
    else:
        raise ValueError

    return abscissa, ordinates


# служебная функция для получения описания МНК
def _get_lsm_description(
        abscissa: list[float], ordinates: list[float]
) -> LSMDescription:
    global event_logger, PRECISION
    n = len(abscissa)
    xy = sum(abscissa[i] * ordinates[i] for i in range(n)) / n
    x_2 = sum(pow(abscissa[i], 2) for i in range(n)) / n
    x = sum(abscissa) / n
    y = sum(ordinates) / n
    a = (xy - x * y) / (x_2 - x ** 2)
    b = y - a * x

    sum_sigma_y = 0
    for i in range(n):
        sum_sigma_y += (ordinates[i] - a * abscissa[i] - b) ** 2
    sigma_y = (1 / (n - 2)) * sum_sigma_y

    sigma_a = pow((sigma_y / (n * (x_2 - x ** 2))), 0.5)

    sigma_b = pow((sigma_y * x_2) / (n * (x_2 - x ** 2)), 0.5)

    return LSMDescription(
        incline=a,
        shift=b,
        incline_error=sigma_a,
        shift_error=sigma_b
    )
