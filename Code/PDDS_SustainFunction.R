
# S1 (first score): any index (sliding window)
# S2 (second score): first observed change at least X months after baseline/S1
# S3 (third score): First observed VALUE (not necessarily change) at least X months after S2

library(pacman)
p_load(dplyr, lubridate)

pdds_sustain = function(data
                        , date = "date"
                        , score = "score"
                        , id = "PATIENT_NUM"
                        , .cleanup = T
                        , .indx = F
                        , window  = dmonths(6)){
  
  `%nin%` = Negate(`%in%`)
  if(any( c(date, score, id) %nin% colnames(data) ) ) stop("Missing date, score, and/or ID column from data frame")
  
  df = data
  colnames(df)[colnames(df) == date] = "date"
  colnames(df)[colnames(df) == score] = "score"
  colnames(df)[colnames(df) == id] = "id"
  
  for(i in 1:nrow(df)){
    
    S1 = S2_index = S3_index = NA 
    # Step 1: Determine the starting score (S1)
    S1 <- df$score[i]
    
    # Step 2: First score with some change after 6 months 
    S2_index <- which((df$score >= S1 + 1 | df$score <= S1 - 1) & df$date > df$date[i] + window & df$id == df$id[i])[1]
    S2 = df$score[S2_index]
    
    # Step 3: Subsequent first score 6 months after S2 date 
    S3_index = which(df$date > df$date[S2_index] + window & df$id == df$id[i])[1]
    S3 = df$score[S3_index]  
    
    # Determine sustainment and direction
    if(is.na(S3)){
      df[i,'Sustainment'] = "Insufficient data"
    }else if(S2-S1>0 & S3-S2>=0){
      df[i,'Sustainment'] = "Progression/Increase"
    }else if(S2-S1 < 0 & S3-S2<=0){
      df[i,'Sustainment'] = "Improvement/Decrease"
    }else{
      df[i,'Sustainment'] = "No sustainment"
    }
    
  }
  
  if(!.cleanup) {
    colnames(df)[colnames(df) == "date"] = date
    colnames(df)[colnames(df) == "score"] = score
    colnames(df)[colnames(df) == "id"] = id
    
    return(df)
  }
  # cleaning up 
  
  df$check = as.integer(df$Sustainment %in% c("Progression/Increase", "Improvement/Decrease"))
  
  df.clean = df %>% group_by(id) %>%
    arrange(id, date) %>% 
    mutate(index = row_number()) %>% 
    arrange(id, desc(check), date) %>% 
    mutate(Sustainment = first(Sustainment)
           , index = first(index)) %>% 
    arrange(id, date) %>% ungroup()
  
  
  if(!.indx) df.clean$index = NULL
  
  df.clean$check = NULL 
  
  
  
  
  colnames(df.clean)[colnames(df.clean) == "date"] = date
  colnames(df.clean)[colnames(df.clean) == "score"] = score
  colnames(df.clean)[colnames(df.clean) == "id"] = id
  
  return(df.clean)
  
  
}


pdds_average_fn = function(data
                        , date = "date"
                        , score = "score"
                        , id = "PATIENT_NUM"
                        , .cleanup = T
                        , .indx = F
                        , window  = dmonths(6)){
  
  `%nin%` = Negate(`%in%`)
  if(any( c(date, score, id) %nin% colnames(data) ) ) stop("Missing date, score, and/or ID column from data frame")
  
  df = data
  colnames(df)[colnames(df) == date] = "date"
  colnames(df)[colnames(df) == score] = "score"
  colnames(df)[colnames(df) == id] = "id"
  
  for(i in 1:nrow(df)){
    rm_index <- which(df$date > df$date[i] & (df$date < df$date[i] + window) & df$id == df$id[i])

    if(length(rm_index)>0) df = df[-c(rm_index),]
  }
  
  if(!.cleanup) {
    colnames(df)[colnames(df) == "date"] = date
    colnames(df)[colnames(df) == "score"] = score
    colnames(df)[colnames(df) == "id"] = id
    
    return(df)
  }
  # cleaning up 
  
  
  
  df.clean = df %>% 
    group_by(id) %>%
    mutate(Average_PDDS = mean(score))

  colnames(df.clean)[colnames(df.clean) == "date"] = date
  colnames(df.clean)[colnames(df.clean) == "score"] = score
  colnames(df.clean)[colnames(df.clean) == "id"] = id
  
  return(df.clean)
  
  
}


