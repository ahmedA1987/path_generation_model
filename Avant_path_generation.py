"""
Name : - Ahmed Abdulsahib
Position: - Research Assistant
Organization: - Tampere University, Beam Research team
Supervisors: Jukka, Reza

Implementation Information
    1) Purpose: -
        This python script represent the model of Embedded Error space generation
        and the purpose is to calculate the Error in distance between the desired position on the
        path where the robot must be and the actual robot position. also the error in Heading
        between the real robot heading and the desired heading at the desired position.
    2) path description:
        the path is calculated using 5th order polynomial Bezier Curve (given 6 control points)
"""
import numpy as np
import math


def binomial_coefficient(bezier_degree, point_number):
    """
    bc = factorial(n)/(factorial(i) * factorial(n - i)) equation from simulink

    :param bezier_degree:
    :param point_number:
    :return:
    """
    coefficient_calculation = math.factorial(bezier_degree) / (math.factorial(point_number) *
                                                               math.factorial(bezier_degree - point_number))
    return coefficient_calculation


def bernstein_polynomial(bezier_degree, point_number, S):
    """
    bp = BinomialCoefficient(n, i) * (1 - t)^(n - i) * t^i; equation from simulink

    :param bezier_degree:
    :param point_number:
    :param S:
    :return:
    """
    coefficient_value = binomial_coefficient(bezier_degree, point_number)
    polynomial_calculation = coefficient_value * pow(1 - S, bezier_degree - point_number) * pow(S, point_number)
    return polynomial_calculation


# PathGen(vh,x,y,psiB,g,drive_cmd,pc104_g_end)
def bezier_curve(control_points, S, desired_point_shift):
    """
    this function calculate the point_number where the robot must be on the path at S
    and calculate the derivative at that point_number also the derivative at point_number slightly different from the
    desired point_number

    :param control_points: type list of arrays: list of the control points
    :param S: type float: value between 0 (start of the path) and 1 (end of the path)
    :param desired_point_shift: type float: the shift of the desired point_number to get second point_number close to
                                the first
    :return:
    math.factorial(5)
    desired_point, derivative_desired_point, derivative_sdp---> P,Pp,P_del
    bc = bc +   BernsteinPolynomial(n, i, t) * P(:,i+1); equation from simulink
    """
    bezier_degree = len(control_points) - 1
    desired_point = np.array([[0], [0]])  # bc
    derivative_desired_point = np.array([[0], [0]])  # bv
    derivative_sdp = np.array([[0], [0]])  # bv_del
    for point_number in range(len(control_points)):
        polynomial_value = bernstein_polynomial(bezier_degree, point_number, S)
        desired_point = desired_point + bezier_degree * polynomial_value * control_points[point_number]
        """
        if i~=n,
        bv = bv + n*BernsteinPolynomial(n-1, i, t) * (P(:,i+2)-P(:,i+1));
        bv_del = bv_del + n*BernsteinPolynomial(n-1, i, t+del) * (P(:,i+2)-P(:,i+1));
        """
        if point_number == bezier_degree - 1:
            polynomial_value = bernstein_polynomial(bezier_degree - 1, point_number, S)
            derivative_desired_point = derivative_desired_point + bezier_degree * polynomial_value * \
                                       (control_points[point_number + 1] - control_points[point_number])

            polynomial_value = bernstein_polynomial(bezier_degree - 1, point_number, S + desired_point_shift)
            derivative_sdp = derivative_sdp + bezier_degree * polynomial_value * \
                                       (control_points[point_number + 1] - control_points[point_number])

    return desired_point, derivative_desired_point, derivative_sdp


def path_generation(heading_velocity, robot_x, robot_y, robot_heading, S, drive_command, S_max=0.99):
    """
    this function will calculate the errors in position and heading between the desired point on the path and
    the robot real position
    :param heading_velocity: type float: speed of the robot
    :param robot_x: type float: robot position on the x axis in world frame
    :param robot_y: type float: robot position on the y axis in world frame
    :param robot_heading: type float: the angle between the x axis of the robot frame and x axis of world frame
    :param S: type float: value between 0 (at begging of the path) and 1 (at end of the path)
    :param drive_command: type vector[floats]:vector has 6 control points and 4 velocities (different velocity for
                                              every 25% of the path) and other unused parameters in this model
    :param S_max: type integer: equals 1, and simply it is the highest value that S can reach
        # Dg, Ds, vr, ye, psiE, Cc,distance_to_go
    :return S_derivative_normalize: value between 0 and 1 that is calculated from normalizing S_derivative
    :return S_derivative: represent the finished distance from  the path measured in meters
    :return reference_velocity: is the desired velocity that is calculated depending on the 4 velocities from the
                            command line ( the equations will be explained in the code below)
    :return heading_error: the differance between the robot angle and the desired angle
    :return y_distance_error: the distance between the robot real location and the desired location on the path on
                            y axis
    :return bezier_curvature: curvature of Bezier Curve
    :return Distance_to_go: the remaining distance
    """

    robot_position = np.array([[robot_x], [robot_y]])  # the robot position as array shape 2x1
    control_constant = 2  # this constant will be multiplied by the saturation in S_derivative equation
    x_error_limit = 1  # the distance differance limit between the robot position and desired position on x axis

    # there are 6 control points from drive command vector will be used to calculate Bezier Curve
    control_point_0 = np.array([[drive_command[0]], [drive_command[1]]])
    control_point_1 = np.array([[drive_command[2]], [drive_command[3]]])
    control_point_2 = np.array([[drive_command[4]], [drive_command[5]]])
    control_point_3 = np.array([[drive_command[6]], [drive_command[7]]])
    control_point_4 = np.array([[drive_command[8]], [drive_command[9]]])
    control_point_5 = np.array([[drive_command[10]], [drive_command[11]]])

    # the reference velocity is set to 0 as initial value
    reference_velocity = 0

    # the 4 velocities from the drive command vector will be used to determined reference velocity
    velocity_0 = drive_command[12]
    velocity_1 = drive_command[13]
    velocity_2 = drive_command[14]
    velocity_3 = drive_command[15]

    # determining the heading speed direction depends on the velocity_0 value as follows
    if velocity_0 < 0:  # meaning the robot is moving backward
        heading_velocity = 0 - heading_velocity

    # desired_point_shift is a small value to calculate the tangent to the Bezier Curve at the
    # desired position for calculating the desired heading
    desired_point_shift = 0.1
    # changing the sign to negative when the desired point is at the end of the curve to calculate the tangent
    if S > 1 - desired_point_shift:
        desired_point_shift = 0 - desired_point_shift

    # desired_point is the position where the robot must be on the path
    # derivative_desired_point is the derivative at the desired point
    # derivative_sdp is the derivative at point slightly different position from the desired point
    desired_point, derivative_desired_point, derivative_sdp = bezier_curve([control_point_0, control_point_1,
                                                                            control_point_2, control_point_3,
                                                                            control_point_0], S, desired_point_shift)









"""
A = [1;2]
B = [3;4]
C = [5;6]
D = [7;8]
E = [9;10]
F = [11;12]
k=[A B C D E F]
s =  size(k, 2)
n = size(k, 2) - 1 

import numpy as np
p0 = np.array([[1], [2]])
p1 = np.array([[3], [4]])
p2 = np.array([[5], [6]])
p3 = np.array([[7], [8]])
p4 = np.array([[9], [10]])
p5 = np.array([[11], [12]])
k = np.concatenate((p0, p1, p2, p3, p4, p5), axis=1) #  2x6
print(k.shape)
"""
