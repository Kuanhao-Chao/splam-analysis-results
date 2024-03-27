            junction_tp = j_pred[np.bitwise_and((j_pred >= THRESHOLD), (j_label == 1))]
            junction_tn = j_pred[np.bitwise_and((j_pred < THRESHOLD), (j_label == 0))]
            junction_fp = j_pred[np.bitwise_and((j_pred >= THRESHOLD), (j_label == 0))]
            junction_fn = j_pred[np.bitwise_and((j_pred < THRESHOLD), (j_label == 1))]

            print("junction_sum : ", len(junction_tp) + len(junction_tn) + len(junction_fp) + len(junction_fn))
            print("junction_tp: ", len(junction_tp))
            print("junction_tn: ", len(junction_tn))
            print("junction_fp: ", len(junction_fp))
            print("junction_fn: ", len(junction_fn))

            print("junction precision: ", len(junction_tp) / (len(junction_tp) + len(junction_fp)))
            print("junction recall   : ", len(junction_tp) / (len(junction_tp) + len(junction_fn)))
            print("junction accuracy : ", (len(junction_tp) + len(junction_tn)) / (len(junction_tp) + len(junction_tn) + len(junction_fp) + len(junction_fn)))