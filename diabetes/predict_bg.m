function bg = predict_bg(M0, S, Ai, e, dt, F, I)
% predict_bg(M0, S, Ai, e, dt, F, I)
% Mi = M0 - S (Ai + e * dt) + (I*F);
% M0 - initial blood glucose level
% S  - sensitivity (BG/unit insulin)
% Ai - amount of insulin that has entered the blood stream
% e  - error in basal (current_basal - true basal)
% dt - time since M0 was measured
% F  - carbohydrates that have entered the blood stream since M0/dt
% I  - BG sensitivity to carbs (amount BG is raised/unit carb)

bg = M0 - S * (Ai + e * dt) + (I*F) ;

