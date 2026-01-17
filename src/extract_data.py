import cv2
import numpy as np
import pandas as pd
import math
import os

# Configuration
VIDEO_PATH = os.path.join(os.path.dirname(__file__), "..", "res", "First_Video_2s.mp4")
OUTPUT_CSV = os.path.join(os.path.dirname(__file__), "..", "res", "video_data.csv")
TARGET_HEIGHT = 800


def resize_frame(frame, target_h):
    """
    Resize a video frame to a specified height while maintaining the aspect ratio.

    Args:
        frame (numpy.ndarray): The input video frame.
        target_h (int): The target height for the resized frame.

    Returns:
        numpy.ndarray: The resized video frame.
    """
    h, w = frame.shape[:2]
    aspect_ratio = w / h
    target_w = int(target_h * aspect_ratio)
    return cv2.resize(frame, (target_w, target_h))


def calculate_angle(p_origin, p_target):
    """
    Calculate the angle between two points relative to the vertical axis.

    Args:
        p_origin (tuple): The origin point (x, y).
        p_target (tuple): The target point (x, y).

    Returns:
        float: The angle in radians.
    """
    dx = p_target[0] - p_origin[0]
    dy = p_origin[1] - p_target[1]  # Inverted y-axis for image coordinates
    return math.atan2(dx, -dy)


def dist(p1, p2):
    """
    Calculate the Euclidean distance between two points.

    Args:
        p1 (tuple): The first point (x, y).
        p2 (tuple): The second point (x, y).

    Returns:
        float: The Euclidean distance between the two points.
    """
    return math.sqrt((p1[0] - p2[0]) ** 2 + (p1[1] - p2[1]) ** 2)


def main():
    # Check if the video file exists
    if not os.path.exists(VIDEO_PATH):
        print(f"ERROR: Video not found: {VIDEO_PATH}")
        exit()

    # Open the video file
    cap = cv2.VideoCapture(VIDEO_PATH)
    fps = cap.get(cv2.CAP_PROP_FPS)
    ret, frame_raw = cap.read()
    if not ret:
        print("Error reading video.")
        exit()

    frame = resize_frame(frame_raw, TARGET_HEIGHT)

    # Select the pivot point
    print("STEP 1: Click on the PIVOT -> Press SPACE/ENTER")
    bbox_pivot = cv2.selectROI("Tracker", frame, False)
    pivot_pos = (int(bbox_pivot[0] + bbox_pivot[2] / 2), int(bbox_pivot[1] + bbox_pivot[3] / 2))

    # Select an orange marker
    print("STEP 2: Select an ORANGE MARKER -> Press SPACE/ENTER")
    bbox_color = cv2.selectROI("Tracker", frame, False)
    x, y, w, h = bbox_color
    roi_color = frame[y:y + h, x:x + w]

    # Calibration: Detect the average color of the selected region
    hsv_roi = cv2.cvtColor(roi_color, cv2.COLOR_BGR2HSV)
    mean_hsv = np.mean(hsv_roi, axis=(0, 1))
    print(f"Detected color: {mean_hsv}")

    # Define the lower and upper bounds for the orange color in HSV
    lower_orange = np.array([max(0, int(mean_hsv[0]) - 20), 50, 50], dtype=np.uint8)
    upper_orange = np.array([min(180, int(mean_hsv[0]) + 20), 255, 255], dtype=np.uint8)

    # Initialize data storage and frame index
    data = []
    frame_idx = 0

    # Memory of last positions
    prev_m1_pos = None
    prev_m2_pos = None

    print("Processing")

    while True:
        ret, frame_raw = cap.read()
        if not ret:
            break

        frame = resize_frame(frame_raw, TARGET_HEIGHT)

        # Create a mask for the orange color
        hsv = cv2.cvtColor(frame, cv2.COLOR_BGR2HSV)
        mask = cv2.inRange(hsv, lower_orange, upper_orange)

        # Clean the mask using erosion and dilation
        mask = cv2.erode(mask, None, iterations=2)
        mask = cv2.dilate(mask, None, iterations=2)

        # Find contours in the mask
        contours, _ = cv2.findContours(mask, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)

        # Keep contours larger than 10 pixels to avoid noise
        valid_contours = [c for c in contours if cv2.contourArea(c) > 10]

        # Sort contours by size (largest to smallest)
        valid_contours.sort(key=cv2.contourArea, reverse=True)

        # Keep the two largest contours (Mass 1 and Mass 2)
        top_contours = valid_contours[:2]

        centers = []
        for c in top_contours:
            M = cv2.moments(c)
            if M["m00"] != 0:
                cx = int(M["m10"] / M["m00"])
                cy = int(M["m01"] / M["m00"])
                centers.append((cx, cy))

        # Assign logic for the two masses
        if len(centers) == 2:
            c1 = centers[0]
            c2 = centers[1]

            m1_pos = None
            m2_pos = None

            if prev_m1_pos is None:
                # First frame: assign based on distance to pivot
                if dist(c1, pivot_pos) < dist(c2, pivot_pos):
                    m1_pos, m2_pos = c1, c2
                else:
                    m1_pos, m2_pos = c2, c1
            else:
                # Subsequent frames: assign based on proximity to previous positions
                if dist(c1, prev_m1_pos) + dist(c2, prev_m2_pos) < dist(c2, prev_m1_pos) + dist(c1, prev_m2_pos):
                    m1_pos, m2_pos = c1, c2
                else:
                    m1_pos, m2_pos = c2, c1

            prev_m1_pos = m1_pos
            prev_m2_pos = m2_pos

            # Calculate angles [rad]
            theta1 = calculate_angle(pivot_pos, m1_pos)
            theta2 = calculate_angle(m1_pos, m2_pos)

            # Calculate length of rods [pixels]
            l1 = dist(pivot_pos, m1_pos)
            l2 = dist(m1_pos, m2_pos)

            time_s = frame_idx / fps
            data.append([time_s, theta1, theta2, l1, l2])

            # Draw lines and circles for visualization
            cv2.line(frame, pivot_pos, m1_pos, (255, 0, 0), 2)
            cv2.line(frame, m1_pos, m2_pos, (0, 0, 255), 2)

            cv2.circle(frame, m1_pos, 8, (0, 255, 0), 2)
            cv2.circle(frame, m2_pos, 8, (0, 255, 0), 2)
        else:
            # If not exactly 2 masses are found, draw the detected ones in red (debug)
            for p in centers:
                cv2.circle(frame, p, 5, (0, 0, 255), -1)

        cv2.imshow("Tracker", frame)
        cv2.imshow("Mask", mask)

        if cv2.waitKey(1) & 0xFF == ord('q'):
            break

        frame_idx += 1

    cap.release()
    cv2.destroyAllWindows()

    # Export the data to a CSV file
    if data:
        df = pd.DataFrame(data, columns=["time_s", "theta1_exp", "theta2_exp", "l1_exp", "l2_exp"])
        df.to_csv(OUTPUT_CSV, index=False, sep=';')
        print(f"Done! {len(data)} points saved in '{OUTPUT_CSV}'.")
    else:
        print("No valid data extracted.")


if __name__ == "__main__":
    main()