import numpy as np
import math

Gamma = 1.4
R = 8.31451

def calculate_theta_degrees(beta_degrees, M1, gamma):

    # beta を度数からラジアンに変換
    beta_rad = np.radians(beta_degrees)

    # sin(beta_rad)がゼロの場合、cot(beta)は未定義になります。
    # ほぼゼロの場合も考慮して、非常に小さい値と比較します。
    if np.isclose(np.sin(beta_rad), 0):
        print(f"エラー: betaが{beta_degrees}度（0または180度の倍数）のため、cot(beta)が未定義です。")
        return None

    cot_beta = np.cos(beta_rad) / np.sin(beta_rad)

    numerator = M1**2 * np.sin(beta_rad)**2 - 1
    denominator = M1**2 * (gamma + np.cos(2 * beta_rad)) + 2

    # 分母がゼロになる可能性をチェック
    # 非常に小さい値（ゼロに近い値）も考慮します。
    if np.isclose(denominator, 0):
        print("エラー: 計算の分母がゼロに非常に近いため、結果が無限大になる可能性があります。")
        return None

    # tan(theta) の値
    tan_theta_val = 2 * cot_beta * (numerator / denominator)

    # theta をラジアンで計算
    theta_rad = np.arctan(tan_theta_val)

    # theta をラジアンから度数に変換
    theta_degrees = np.degrees(theta_rad)

    return theta_degrees

def calculate_Mnext(gamma, M1, beta_degrees, theta_degrees):

    # 度数をラジアンに変換
    beta_radians = math.radians(beta_degrees)
    theta_radians = math.radians(theta_degrees)

    # sin(β - θ) の計算
    sin_beta_minus_theta = math.sin(beta_radians - theta_radians)

    # ゼロ除算のチェック
    if abs(sin_beta_minus_theta) < 1e-9:  # 非常に小さい値もゼロとみなす
        raise ValueError("sin(beta - theta)がゼロに非常に近いため、M2は計算できません。これは通常、物理的に不可能な状態を示します。")

    # 式の左側の部分
    term_left = 1 / sin_beta_minus_theta

    # 平方根の中の分子
    numerator_sqrt = 1 + ((gamma - 1) / 2) * M1**2 * math.sin(beta_radians)**2

    # 平方根の中の分母
    denominator_sqrt = gamma * M1**2 * math.sin(beta_radians)**2 - (gamma - 1) / 2




    # 式の右側の平方根の部分
    term_right = math.sqrt(numerator_sqrt / denominator_sqrt)

    # M2 の最終計算
    M2 = term_left * term_right
    return M2

def calculate_p02_p01(gamma, M1, beta_degrees):

    beta_radians = math.radians(beta_degrees)  # 度数をラジアンに変換
    sin_squared_beta = math.sin(beta_radians)**2

    term1_base = 1 + (2 * gamma / (gamma + 1)) * (M1**2 * sin_squared_beta - 1)
    term1_exponent = -1 / (gamma - 1)
    term1 = term1_base**term1_exponent

    term2_numerator = 2 + (gamma - 1) * M1**2 * sin_squared_beta
    term2_denominator = (gamma + 1) * M1**2 * sin_squared_beta
    term2_base = term2_numerator / term2_denominator
    term2_exponent = -gamma / (gamma - 1)
    term2 = term2_base**term2_exponent

    p02_p01 = term1 * term2
    return p02_p01

def D_entropy(PRatio):
    answer = (-1) * math.log(PRatio,math.e)
    return answer

def case1():

    Mn = [0.0, 4.0, 0.0, 0.0, 0.0, 0.0]
    Beta = [0.0, 20, 30, 50, 90, 0.0]

    Theta_Result_Deg = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    PRatio_Result = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    D_entropy_result = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    Addition_entropy_result = 0.0
    Addition_PRatio_result = 1

    for i in range(1,5):

        #theata算出/代入
        Theta_Result_Deg[i+1] = calculate_theta_degrees(Beta[i], Mn[i], Gamma)

        #M.next計算
        Mn[i+1] = calculate_Mnext(Gamma, Mn[i], Beta[i], calculate_theta_degrees(Beta[i], Mn[i], Gamma))

        #p.ratio計算
        PRatio_Result[i+1] = calculate_p02_p01(Gamma, Mn[i], Beta[i])

        #D.entropyn計算
        D_entropy_result[i+1] = D_entropy(PRatio_Result[i+1])


        #Addition_d_entropy計算
        Addition_entropy_result = Addition_entropy_result + D_entropy_result[i+1]

        #Addition_PRatio_resul計算
        Addition_PRatio_result = Addition_PRatio_result * PRatio_Result[i+1]

    print('Case1：斜め衝撃波（M.final, P.Ratio.final, Entropy.result）')
    print(Mn[5],', ',Addition_PRatio_result,', ',Addition_entropy_result)

def case2():

    M1 = 4.0
    Beta = 90

    PRatio_Result_C2 = 0.0
    D_entropy_result = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    Addition_entropy_result_C2 = 0.0
    Addition_PRatio_result_C2 = 1
    
    #theata算出/代入
    Theta_Result_Deg_C2 = calculate_theta_degrees(Beta, M1, Gamma)

    #M.next計算
    M2 = calculate_Mnext(Gamma, M1, Beta, calculate_theta_degrees(Beta, M1, Gamma))

    #p.ratio計算
    PRatio_Result_C2 = calculate_p02_p01(Gamma, M1, Beta)

    #D.entropyn計算
    D_entropy_result_C2 = D_entropy(PRatio_Result_C2)

    #Addition_d_entropy計算
    Addition_entropy_result_C2 = Addition_entropy_result_C2 + D_entropy_result_C2

    #Addition_PRatio_resul計算
    Addition_PRatio_result_C2 = Addition_PRatio_result_C2 * PRatio_Result_C2

    print('Case2：垂直衝撃波（M.final, P.Ratio.final, Entropy.result）')
    print(M2,', ',Addition_PRatio_result_C2,', ',Addition_entropy_result_C2)


case1()
case2()